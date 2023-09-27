function [fbstats] = vbhmm_fb(data, hmm_varpar, opt)
% vbhmm_fb - run forward-backward algorithm
%
% this is an internal function called by vbhmm_em
%
% [fbstats] = vbhmm_fb(data, hmm_varpar, opt)
%
% opt.savexi = 1: also save xi for each t (default=0)
% opt.usegroups = 1: use groups mapping (default=0)
% opt.useMEX = 1: use MEX implementation if available (default=1)
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-01-13
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong

% 2017-02-13: ABC: code cleanup (rename some variables). Speed optimization.
% 2017-06-09: ABC: code cleanup (rename some variables).
% 2017-06-09: ABC: added option to calculate derivatives
% 2017-06-14: ABC: code cleanup: epsilon & logATilde changed to row format
% 2017-07-31: ABC: removed calculate derivative options, since it's not used anymore.
% 2017-08-01: ABC: speed-up using MEX implementation
% 2017-08-03: ABC: clean-up code
% 2017-08-25: ABC: added scaling for numerical precision problems
% 2018-03-10: ABC: support groups

% hmm parameters: v, W, epsilon, alpha, m, beta
v       = hmm_varpar.v;
W       = hmm_varpar.W;
epsilon = hmm_varpar.epsilon;
alpha   = hmm_varpar.alpha;
m       = hmm_varpar.m;
beta    = hmm_varpar.beta;

% options: group_map, usegroups
usegroups = opt.usegroups;
if usegroups
  group_inds = opt.group_inds;
  group_map = opt.group_map;
  numgroups = opt.numgroups;
end

% other options: save xi
if isfield(opt, 'savexi')
  savexi = opt.savexi;
else
  savexi = 0;
end
if ~isfield(opt, 'useMEX')
  opt.useMEX = 1;
end

% some constants: K, dim, N, maxT, lengthT, const, const_denominator
K   = size(m, 2);
dim = size(m, 1);
N   = length(data);
maxT    = max(cellfun('prodofsize', data)/dim);
lengthT = maxT;
const = dim*log(2);
const_denominator = (dim*log(2*pi))/2;

% pre-calculate constants
logLambdaTilde = zeros(1,K);
for k = 1:K
  t1 = psi(0, 0.5*repmat(v(k)+1,dim,1) - 0.5*[1:dim]');
  logLambdaTilde(k) = sum(t1) + const  + log(det(W(:,:,k)));
end

if ~usegroups
  psiEpsilonHat = zeros(1,K);
  logATilde     = zeros(K,K);
  for k = 1:K
    %Bishop (10.66)
    psiEpsilonHat(k) = psi(0,sum(epsilon(k,:)));
    logATilde(k,:) = psi(0,epsilon(k,:)) - psiEpsilonHat(k);  % A(i,j) = p(j->i) [column]
  end
  psiAlphaHat = psi(0,sum(alpha));
  logPiTilde = psi(0,alpha) - psiAlphaHat;
  
else
  psiEpsilonHat = {};   logATilde = {};
  psiAlphaHat = {};     logPiTilde = {};
  for g = 1:numgroups
    for k = 1:K
      %Bishop (10.66)
      psiEpsilonHat{g}(k) = psi(0,sum(epsilon{g}(k,:)));
      logATilde{g}(k,:) = psi(0,epsilon{g}(k,:)) - psiEpsilonHat{g}(k);  % A(i,j) = p(j->i) [column]
    end
    psiAlphaHat{g} = psi(0,sum(alpha{g}));
    logPiTilde{g} = psi(0,alpha{g}) - psiAlphaHat{g};
  end
end


%% check if MEX can be used
% MEX does not support savexi
canuseMEX = 0;
if (opt.useMEX)
  if (exist('vbhmm_fb_mex')==3)
    canuseMEX = 1;
    if (savexi)
      canuseMEX = 0;
      warning('vbhmm_fb: savexi not supported in MEX implementation');
    end
    %if (usegroups)
    %  canuseMEX = 0;
    %  warning('vbhmm_fb: usegroups not supported in MEX implementation');
    %end
  else
    canuseMEX = 0;
    warning('vbhmm_fb: MEX implementation not found');
  end
end

if (opt.useMEX) && (canuseMEX)
  %% use MEX implementation
  
  if ~usegroups  
    % setup some variables first
    t_pz1     = exp(logPiTilde'); % priors p(z1)
    t_tpztzt1 = exp(logATilde); % transitions p(zt|z(t-1)) [row format]
    
    % Inputs
    %   data = {Nx1} - data{i} = [dim x T]
    %   K, N, dim, maxT = scalars [1x1]
    %   m = [dim x K]
    %   W = [dim x dim x K]
    %   v = [K x 1]
    %   beta = [K x 1]
    %   logLambdaTilde = [1 x K]
    %   const_denominator = scalar [1 x 1]
    %   t_pz1 = [1 x K]
    %   t_tpz1zt1 = [K x K]
    %
    % Outputs:
    %   logrho_Saved = [K x N x maxT]
    %   gamma_all    = [K x N x maxT]
    %   xi_sum       = [K x K x N]
    %   phi_norm     = [1 x N]
    %   xi_Saved (unsupported, since currently unused elsewhere)
    
    % call MEX function
    [logrho_Saved, gamma_all, xi_sum, phi_norm] = ...
      vbhmm_fb_mex(data, K, N, dim, maxT, m, W, v, beta, logLambdaTilde, const_denominator, t_pz1, t_tpztzt1);
    
  else
    % use groups
    
    % storage
    logrho_Saved = zeros(K, N, lengthT);
    phi_norm = zeros(1,N);
    gamma_all = zeros(K,N,maxT);
    xi_sum = zeros(K,K,N);
  
    for g=1:numgroups
      % run forward-backward for data in each group
      g_inds = group_inds{g};
      
      % get params for this group
      t_pz1      = exp(logPiTilde{g}');  % priors p(z1)
      t_tpztzt1  = exp(logATilde{g});    % transitions p(zt|z(t-1)) [row format]
      
      % select data for this group
      g_data = data(g_inds);
      g_N = length(g_data);
    
      % run MEX file
      [g_logrho_Saved, g_gamma_all, g_xi_sum, g_phi_norm] = ...
        vbhmm_fb_mex(g_data, K, g_N, dim, maxT, m, W, v, beta, logLambdaTilde, const_denominator, t_pz1, t_tpztzt1);
    
      % store results
      logrho_Saved(:,g_inds,:) = g_logrho_Saved;
      gamma_all(:,g_inds,:)    = g_gamma_all;
      phi_norm(g_inds)         = g_phi_norm;
      xi_sum(:,:,g_inds)       = g_xi_sum;
    end
    
    % check w/ non-MEX version    
    %tmpopt = opt;
    %tmpopt.useMEX = 0;
    %tmpfbstats = vbhmm_fb(data, hmm_varpar, tmpopt);
    %
    %e1 = sum(abs(tmpfbstats.logrho_Saved(:) - logrho_Saved(:)));
    %e2 = sum(abs(tmpfbstats.gamma_all(:) - gamma_all(:)));
    %e3 = sum(abs(tmpfbstats.phi_norm(:) - phi_norm(:)));
    %e4 = sum(abs(tmpfbstats.xi_sum(:) - xi_sum(:)));
    %
    %tmp = [e1, e2, e3, e4];
    %if any(tmp > 1e-9)
    %  keyboard
    %end    
  end
    
else
  if (opt.useMEX)
    % message if we wanted to use MEX but couldn't
    warning('vbhmm_fb: using MATLAB implementation. This will be slower');
  end
  
  %% MATLAB implementation
  
  % initialize storage
  logrho_Saved = zeros(K, N, lengthT);
  phi_norm = zeros(1,N);
  gamma_all = zeros(K,N,maxT);
  xi_sum = zeros(K,K,N);
  if savexi
    xi_Saved = cell(1,N);
    for n=1:N
      xi_Saved{n} = zeros(K,K,size(data{n},1));
    end
  end
  
  % precompute pz1 and tpztzt1
  if ~usegroups
    t_pz1     = exp(logPiTilde'); % priors p(z1)
    t_tpztzt1 = exp(logATilde);   % transitions p(zt|z(t-1)) [row format]
  else
    for g=1:numgroups   
      gt_pz1{g}         = exp(logPiTilde{g}');  % priors p(z1)
      gt_tpztzt1{g}     = exp(logATilde{g});    % transitions p(zt|z(t-1)) [row format]
    end    
  end
  
  %% run FB on each sequence
  for n = 1:N
    tdata = data{n}; tdata = tdata';
    tT = size(tdata,2);
    %delta = []; logrho = [];
    
    delta = zeros(K,tT);
    
    for k = 1:K
      %OLD slow code Bishop (10.64)
      %for t = 1:tT
      %  diff = tdata(:,t) - m(:,k);
      %  delta(k,t) = dim/beta(k) + v(k)*diff'*W(:,:,k)*diff;
      %end
      %delta_old = delta;
      
      % ABC: 2016-04-21 - fast code
      diff = bsxfun(@minus, tdata, m(:,k));
      mterm = sum((W(:,:,k)*diff).*diff,1);
      delta(k,:) = dim/beta(k) + v(k) * mterm;
      
    end
    % OLD slow code - Bishop (10.46)
    %for k = 1:K
    %  for t = 1:tT
    %    logrho(k,t) = 0.5*logLambdaTilde(k) - 0.5*delta(k,t) - const_denominator;
    %  end
    %end
    %logrho_old = logrho;
    
    % ABC: 2016-04-21 - fast code
    logrho = bsxfun(@minus, 0.5*logLambdaTilde(:), 0.5*delta) - const_denominator;
    
    logrho_Saved(:,n,1:tT) = logrho;
    
    %% forward_backward %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %gamma = zeros(K,tT);
    %sumxi = zeros(K,K);  % [row]
    
    fb_logrho = logrho';
    %if ~usegroups
    %  fb_logPiTilde = logPiTilde';
    %  fb_logATilde = logATilde;
    %else
    %  fb_logPiTilde = logPiTilde{group_map(n)}';
    %  fb_logATilde = logATilde{group_map(n)};
    %end
    %t_alpha = [];t_beta = []; t_c = [];
    t_alpha = zeros(tT, K); % alpha hat
    t_Delta = zeros(tT, K); % Delta (unnormalized version)
    t_beta  = zeros(tT, K);
    t_Eta   = zeros(tT, K); % Eta (unnormalized version)
    t_c     = zeros(1,tT);
    
    %t_pz1 = exp(fb_logPiTilde);    % priors p(z1)
    %t_tpztzt1 = exp(fb_logATilde); % transitions p(zt|z(t-1)) [row format]
    if usegroups
      t_pz1 = gt_pz1{group_map(n)};
      t_tpztzt1 = gt_tpztzt1{group_map(n)};
    end
    
    %t_pxtzt = exp(fb_logrho);      % emissions p(xt|zt)
    
    % emission p(xt|zt)/max(p(xt|zt)) (more stable)
    max_fb_logrho = max(fb_logrho, [], 2);
    t_pxtzt = exp(bsxfun(@minus, fb_logrho, max_fb_logrho));
    
    t_x = tdata';
    t_T = size(t_x,1);
    
    t_gamma = zeros(K,t_T);
    t_sumxi = zeros(K,K);
      
    if t_T >= 1
      %forward
      %t_gamma = [];
      
      %before normalizing
      t_Delta(1,:) = t_pz1.*t_pxtzt(1,:);
      
      % 2016-04-29 ABC: rescale for numerical stability (otherwise values get too small)
      t_c(1) = sum(t_Delta(1,:));
      
      % normalize
      t_alpha(1,:) = t_Delta(1,:) / t_c(1);
      
      if t_T > 1
        for i=2:t_T
          % before normalizing
          t_Delta(i,:) = t_alpha(i-1,:)*t_tpztzt1.*t_pxtzt(i,:);
          
          % 2016-04-29 ABC: rescale for numerical stability
          t_c(i) = sum(t_Delta(i,:));
          
          % normalize
          t_alpha(i,:) = t_Delta(i,:) / t_c(i);
        end
      end
      
      %% backward %%%%%%%%%%%%%%%%%%%%%%
      % BUG FIX 2017-02-09: ones(1,K)./K;
      % doesn't affect HMM results, since we had normalized anyways.
      t_beta(t_T,:) = ones(1,K);
      
      t_gamma(:,t_T) = (t_alpha(t_T,:).*t_beta(t_T,:))';
      
      if t_T > 1
        for i=(t_T-1):-1:1
          % before normalizing
          bpi = (t_beta(i+1,:).*t_pxtzt(i+1,:));
          t_Eta(i,:) = bpi*t_tpztzt1';
          
          % 2016-04-29 ABC: rescale
          t_beta(i,:) = t_Eta(i,:)/t_c(i+1);
          
          % responsibility
          t_gamma(:,i) = (t_alpha(i,:).*t_beta(i,:))';
          
          %t_sumxi = t_sumxi + (t_logATilde.*(t_alpha(i,:)'*bpi));
          
          % 2016-04-29 ABC: rescale xi
          tmp_xi = (t_tpztzt1.*(t_alpha(i,:)'*bpi)) / t_c(i+1);
          
          % 2016-04-29 ABC BUG FIX: normalize xi matrix to sum to 1
          % (it's a joint probability matrix)
          % 2017-02-09 - not necessary anymore after fixing beta initialization bug.
          %tmp_xi = tmp_xi / sum(tmp_xi(:));
          
          % accumulate
          t_sumxi = t_sumxi + tmp_xi;
          
          if savexi
            xi_Saved{n}(:,:,i) = tmp_xi;
          end
          
        end
      end
      %for i = 1:size(t_gamma,2)
      %  gamma(:,i) = t_gamma(:,i);
      %end
      %gamma(:,1:tT) = t_gamma;      
      %sumxi = t_sumxi;
      
    end
    gamma_all(:,n,1:tT) = t_gamma(:,1:tT);
    xi_sum(:,:,n) = t_sumxi;
    
    % from scaling constants
    %phi_norm(n) = sum(log(t_c));

    % from scaling constants (more stable)
    phi_norm(n) = sum(log(t_c)) + sum(max_fb_logrho);

  end
end

% output
fbstats.logrho_Saved = logrho_Saved;
fbstats.gamma_all    = gamma_all;
fbstats.xi_sum       = xi_sum;
fbstats.phi_norm     = phi_norm;  %Phi (normalization constant)
fbstats.logLambdaTilde = logLambdaTilde;
fbstats.logPiTilde     = logPiTilde;
fbstats.logATilde      = logATilde;

if savexi
  fbstats.xi_Saved = xi_Saved;
end
