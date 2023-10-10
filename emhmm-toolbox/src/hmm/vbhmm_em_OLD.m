function [hmm,L] = vbhmm_em(data,K,ini)
% vbhmm_em - run EM for vbhmm (internal function)
% Use vbhmm_learn instead. For options, see vbhmm_learn.
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-01-13
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong

% see vbhmm_learn for changes

% ini = vbopt

VERBOSE_MODE = ini.verbose;

% get length of each chain
trial = length(data);
datalen = [];
for i = 1:length(data)
    d = data{i};
    datalen(i,1) = size(d,1);
end
lengthT = max(datalen); % find the longest chain
totalT  = sum(datalen);

% clip hyps at extreme values for more stable computations
[ini, hyp_clipped] = vbhmm_clip_hyps(ini);


%initialize the parameters
mix_t = vbhmm_init(data,K,ini); %initialize the parameters

mix = mix_t;
dim      = mix.dim; % dimension of the data
K        = mix.K;   % no. of hidden states
N        = trial;   % no. of chains
maxT     = lengthT; % the longest chain

alpha0   = mix.alpha0;   % hyper-parameter for the priors
epsilon0 = mix.epsilon0; % hyper-parameter for the transitions
m0       = mix.m0;       % hyper-parameter for the mean
beta0    = mix.beta0;    % hyper-parameter for beta (Gamma)
v0       = mix.v0;       % hyper-parameter for v (Inverse-Wishart)
W0inv    = mix.W0inv;    % hyper-parameter for Inverse-Wishart
W0mode   = mix.W0mode;   % W0 is iid or diag?
alpha    = mix.alpha;    % priors

epsilon  = mix.epsilon;  % transitions [row format, k-th row is conditioned on z_{t-1}=k]
beta     = mix.beta;     % beta (Gamma)
v        = mix.v;        % v (Inverse-Wishart)
m        = mix.m;        % mean
W        = mix.W;        % Inverse-Wishart
C        = mix.C;        % covariance
const    = mix.const;    % constants
const_denominator = mix.const_denominator; %constants
maxIter = ini.maxIter;   % maximum iterations allowed
minDiff = ini.minDiff;   % termination criterion

L = -realmax; %log-likelihood
lastL = -realmax; %log-likelihood

% setup groups
if ~isempty(ini.groups)
  usegroups = 1;
  group_ids = unique(ini.groups);  % unique group ids
  numgroups = length(group_ids);
  for g=1:numgroups
    group_inds{g} = find(ini.groups == group_ids(g)); % indices for group members
  end
  % sanitized group membership (1 to G)
  group_map = zeros(1,length(ini.groups));
  for g=1:numgroups
    group_map(group_inds{g}) = g;
  end
  
  % reshape alpha, epsilon into cell
  % also Nk1 and M are cells
  if ~iscell(epsilon)
    tmp = epsilon;
    tmpa = alpha;
    epsilon = {};
    alpha = {};
    for g = 1:numgroups
      epsilon{g} = tmp;
      alpha{g} = tmpa;
    end
  end
  
else
  usegroups = 0;
end

% calculate derivative?
if ini.calc_LLderiv
  do_deriv = 1;
else
  do_deriv = 0;
end 
  
% make the data block
data_block = cat(1, data{:});

% make indexes into full block
index_block = zeros(1,size(data_block,1));
ind = 0;
for n=1:N
  tT = datalen(n);
  index_block(ind+(1:tT)) = (n-1)*maxT+(1:tT);
  ind = ind+tT;
end

for iter = 1:maxIter
          
    %% E step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% 2016-08-14 - call function for FB algorithm      
    % setup HMM
    fbhmm_varpar.v = v;
    fbhmm_varpar.W = W;
    fbhmm_varpar.epsilon = epsilon;
    fbhmm_varpar.alpha = alpha;
    fbhmm_varpar.m = m;
    fbhmm_varpar.beta = beta;
    fbopt.usegroups = usegroups;
    if usegroups
      fbopt.group_inds = group_inds;
      fbopt.group_map = group_map;
      fbopt.numgroups = numgroups;
    end
    fbopt.useMEX = ini.useMEX;
    
    % call FB algorithm
    [fbstats] = vbhmm_fb(data, fbhmm_varpar, fbopt);
    
    % get statistics and constants
    logrho_Saved  = fbstats.logrho_Saved;
    t_gamma_Saved = fbstats.gamma_all;
    xi_sum = fbstats.xi_sum;
    phi_norm  = fbstats.phi_norm;
    logLambdaTilde = fbstats.logLambdaTilde;
    logPiTilde = fbstats.logPiTilde;
    logATilde = fbstats.logATilde;

    
    %% sum up the responsibilities    
    % ABC: 2016-04-21 - fast code
    % 2017-02-09 - not necessary after minor bug fix in vbhmm_fb
    % gamma is already normalized in vbhmm_fb
    
    %scale = sum(gamma_sum, 1);
    %if any(abs(scale-1)>1e-6)
    %  keyboard
    %end    
    %scale(scale==0) = 1;
    %t_gamma_Saved = bsxfun(@rdivide, gamma_sum, scale);  % marginal responsibilities
    
    % sum over sequences (N)
    t_Nk1 = reshape(sum(t_gamma_Saved, 2), [K, maxT]);
    
    % for updating priors
    if ~usegroups
      Nk1 = t_Nk1(:,1);
      Nk1 = Nk1 + 1e-50;
    else
      for g=1:numgroups
        Nk1{g} = sum(t_gamma_Saved(:,group_inds{g},1), 2);
      end
    end
    
    % ROI counts for updating beta and v
    Nk = sum(t_Nk1,2); 
    Nk = Nk + 1e-50;

    % 2016-04-29: ABC NEW fast code
    % transition counts
    % M(i,j) = soft counts i->j [row]
    if ~usegroups
      M = sum(xi_sum, 3);      
    else
      for g=1:numgroups
        M{g} = sum(xi_sum(:,:,group_inds{g}),3);
      end
    end
        
    % ROI mean fixation
    if 0
      % ABC 2016-04-21 - fast code
      t_xbar = zeros(K,N,maxT,dim);
      for n=1:N
        x = data{n}; tT = size(x,1);
        t_xbar(:,n,1:tT,:) = bsxfun(@times, ...
          reshape(t_gamma_Saved(:,n,1:tT), [K,1,tT,1]), ...
          reshape(x, [1 1 tT dim]));
      end
      
      % sum over time
      % ABC 2016-04-21 - fast code
      t1_xbar = sum(t_xbar,3);
      
      % sum over sequences
      % ABC 2016-04-21 - fast code
      t2_xbar = sum(t1_xbar,2);
      t2_xbar = reshape(t2_xbar, [K, dim]);
      
      % normalize
      % ABC 2016-04-21 - fast code
      xbar = bsxfun(@rdivide, t2_xbar, Nk);
    end
    
    % ABC 2017-08-03 - even faster code
    % make gamma block
    tmp_gamma = reshape(permute(t_gamma_Saved, [1 3 2]), [K, maxT*N]);
    gamma_block = tmp_gamma(:,index_block);
    
    % ROI mean fixation
    xbar = zeros(K,dim);
    for k=1:K
      xbar(k,:) = sum(bsxfun(@times, data_block, gamma_block(k,:)'),1) / Nk(k);
    end
    
    
    % scatter matrix
    if 0
      % ABC 2016-04-21: fast code
      t1_S = zeros(dim, dim, K);
      for n=1:N
        x = data{n}';
        tT = size(x,2);
        
        for k=1:K
          d1 = bsxfun(@minus, x, xbar(k,:)');
          d2 = bsxfun(@times, reshape(t_gamma_Saved(k,n,1:tT), [1,tT]), d1);
          t1_S(:,:,k) = t1_S(:,:,k) + d1*d2';
        end
      end
      t1_S = bsxfun(@rdivide, t1_S, reshape(Nk, [1 1 K]));
    end
    
    % scatter matrix
    % ABC 2017-08-02: even faster code
    t1_S = zeros(dim, dim, K);
    for k=1:K
      d1 = bsxfun(@minus, data_block, xbar(k,:));
      d2 = bsxfun(@times, d1, gamma_block(k,:)');
      t1_S(:,:,k) = d2'*d1 / Nk(k);
    end
    
    %% calculate lower bound
    
    % save previous
    if iter > 1
      lastL = L;
    end
    
    % call LB function
    lbhmm_stats.dim      = dim;
    lbhmm_stats.K        = K;
    lbhmm_stats.N        = N;
    lbhmm_stats.alpha0   = alpha0;
    lbhmm_stats.epsilon0 = epsilon0;
    lbhmm_stats.m0       = m0;
    lbhmm_stats.beta0    = beta0;
    lbhmm_stats.v0       = v0;
    lbhmm_stats.W0inv    = W0inv;
    lbhmm_stats.W0mode   = W0mode;
    lbhmm_stats.usegroups = usegroups;
    if usegroups
      lbhmm_stats.group_inds = group_inds;
    end
    lbhmm_stats.t1_S     = t1_S;
    lbhmm_stats.xbar     = xbar;
    lbhmm_stats.Nk       = Nk;
    lbhmm_stats.M        = M;
    lbhmm_stats.clipped  = hyp_clipped;
    L = vbhmm_em_lb(lbhmm_stats, fbhmm_varpar, fbstats);
            
    do_break = 0;
    unstable = 0;
    
    if iter == 1
      % show the first iteration too
      if (VERBOSE_MODE >= 3)
        fprintf('iter %d: L=%g; dL=x\n', 1, L);
      end
    end
    
    if iter > 1
      likIncr = abs((L-lastL)/lastL);
      if (VERBOSE_MODE >= 3)        
        fprintf('iter %d: L=%g; dL=%g', iter, L, likIncr);
        if (L-lastL < 0)
          fprintf('[LL decreased: %g]!!!', L-lastL);
          %keyboard
        end
        fprintf('\n');
      else 
        if (L-lastL < 0)
          if (VERBOSE_MODE >= 2)
            fprintf('[LL decreased]');
          end
        end
      end
      if (likIncr <= minDiff)
        do_break = 1;
      end
    end
    if (iter == maxIter)
      warning('max iterations reached');
      do_break = 1;
    end
    
    % if L is NaN, then usually it means an impossible model L=-inf.
    % NaN occurs because of a divide by 0 in the FB algorithm.
    if any(isnan(L))
      if (VERBOSE_MODE >=1)
        warning('NaN found, unstable model');
      end
      do_break = 1;
      unstable = 1;
      
      % 2018-02-01: stability fix when optimizing hyps (when it tries weird parameter combos)
      L = -inf;
      
      % save for debugging
      if (VERBOSE_MODE >= 3)
        foo=tempname('.');
        warning('degenerate models detected: saving to %s', foo);
        save([foo '.mat'])
      end
    end

    % calculate derivatives if breaking (do this before the last M-step)
    if (do_break) && (do_deriv)
      lbhmm_stats.do_deriv = 1;
      [Lnew, dL] = vbhmm_em_lb(lbhmm_stats, fbhmm_varpar, fbstats);     

      % invalidate the gradient
      if (unstable)
        tmp = fieldnames(dL);
        for i=1:length(tmp)
          dL.(tmp{i}) = nan*dL.(tmp{i});
        end
      end
    end

    % break early, since unstable model
    if (do_break) && unstable
      break;
    end
      
    
    %% M step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % prior & transition parameters
    if ~usegroups      
      alpha = alpha0 + Nk1;
      epsilon = epsilon0 + M;
    else
      for g=1:numgroups
        alpha{g} = alpha0 + Nk1{g};
        epsilon{g} = epsilon0 + M{g};
      end
    end
    
    % update Gaussians
    if ~ini.fix_clusters
      % mean
      beta = beta0 + Nk;
      v = v0 + Nk + 1;
      for k = 1:K
        m(:,k) = (beta0.*m0 + Nk(k).*xbar(k,:)')./beta(k);
      end
      
      % wishart
      for k = 1:K
        if isempty(ini.fix_cov)
          mult1 = beta0.*Nk(k)/(beta0 + Nk(k));
          diff3 = xbar(k,:) - m0';
          diff3 = diff3';
          W(:,:,k) = inv(W0inv + Nk(k)*t1_S(:,:,k) + mult1*diff3*diff3');
                   
          % 2018-04-16: symmetrize for numerical stability
          W(:,:,k) = (W(:,:,k)+W(:,:,k)')/2;
          
        else
          % 2017-01-21 - fix the covariance matrix (for Antoine)
          % set Wishart W to the appropriate value to get the covariance matrix.
          W(:,:,k) = inv(ini.fix_cov)/(v(k)-dim-1);
        end
      end
    end
    
    % covariance
    for k = 1:K
      % precision ~ Wishart(W, v)
      % hence, covariance ~ inverse Wishart(W^{-1}, v)
      % the mean of inverse Wishart is:
      if (v(k) > dim+1)
        C(:,:,k) = inv(W(:,:,k))/(v(k)-dim-1);
      else
        % 2018-06-15: v0.74 - when mean doesn't exist, use the inverse of the mean of the Wishart
        % this usually only happens when there are no fixations in an ROI
        C(:,:,k) = inv(W(:,:,k)) / v(k);
      end
      
      % 2018-04-16: symmetrize for numerical stability
      C(:,:,k) = (C(:,:,k)+C(:,:,k)')/2;
    end
    
    % break out of E & M loop
    if (do_break)
      break;
    end
end

if (VERBOSE_MODE >= 2)
  fprintf('EM numiter=%d: L=%g; dL=%g\n', iter, L, likIncr);
end

if (VERBOSE_MODE >= 3)
  fprintf('-\n');
end

%generate the output model
% NOTE: if adding a new field, remember to modify vbhmm_permute.
hmm = {};
if ~usegroups
  prior_s = sum(alpha);
  prior = alpha ./ prior_s;
  hmm.prior = prior;
  trans_t = epsilon;  % [row]
  for k = 1:K
      scale = sum(trans_t(k,:));
      if scale == 0 scale = 1; end
      trans_t(k,:) = trans_t(k,:)./repmat(scale,1,K);
  end
  hmm.trans = trans_t;
else
  for g=1:numgroups
    prior_s = sum(alpha{g});
    prior = alpha{g} ./ prior_s;
    hmm.prior{g} = prior;
    
    trans_t = epsilon{g};  % [row]
    for k = 1:K
      scale = sum(trans_t(k,:));
      if scale == 0 scale = 1; end
      trans_t(k,:) = trans_t(k,:)./repmat(scale,1,K);
    end
    hmm.trans{g} = trans_t;
  end
end

hmm.pdf = {};
for k = 1:K
  % Antoine mod - fix the covariance after learning
  %if ~isempty(ini.fix_cov)
  %  %C(:,:,k)=diag(diag(C(:,:,k)));%0s on 12 and 21
  %  C(:,:,k)=ini.fix_cov; %[ini.do_constrain_var_size 0; 0 ini.do_constrain_var_size];
  %end
  
  hmm.pdf{k,1}.mean = m(:,k)';
  hmm.pdf{k,1}.cov = C(:,:,k);
end
hmm.LL = L;
hmm.gamma = cell(1,N);
for n=1:N
  hmm.gamma{n} = reshape(t_gamma_Saved(:,n,1:datalen(n)), [K datalen(n)]);
end
hmm.M = M;     % transition counts
hmm.N1 = Nk1;  % prior counts
hmm.N  = Nk;   % cluster sizes

% save group info
if usegroups
  for g=1:numgroups
    ggamma = hmm.gamma(group_inds{g});
    hmm.Ng{g} = sum(cat(2, ggamma{:}), 2); % cluster size by group
  end
  hmm.group_map = group_map;
  hmm.group_ids = group_ids;
  hmm.group_inds = group_inds;
end

% save variational parameters
hmm.varpar.epsilon = epsilon;
hmm.varpar.alpha = alpha;
hmm.varpar.beta = beta;
hmm.varpar.v = v;
hmm.varpar.m = m;
hmm.varpar.W = W;

% save derivatives
if (do_deriv)
  hmm.dLL = dL;
end
  