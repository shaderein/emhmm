function h3m_r = hem_h3m_c_step(h3m_r,h3m_b,mopt) 
% h3m_r = initialization for the reduced mixture output
% h3m_b = base mixture input
% mopt  = options 
%       .K     = number of mixtures in the reduced model
%       .Nv     = number of virtual samples (this is multiplied by Kb)
%       .tau   = length of virtual sequences
%       .termmode  = how to terminate EM
%       .termvalue = when to terminate the EM
% h3m.K
% h3m.hmm = {hmm1 hmm2 ... hmmK}
% h3m.omega = [omega1 omega2 ... omegaK]
% hmm1.A hmm1.emit hmm1.prior
%
% ---
% H3M Toolbox 
% Copyright 2012, Emanuele Coviello, CALab, UCSD

% 2016-12-09: ABC - added output for emission population size
% 2017-04-25: ABC - small bug fix when generate component is handled.
% 2017-06-06: ABC - put M-step into a separate function
% 2017-08-07: ABC - caching to make E-step faster
% 2017-08-08: ABC - added MEX implementation for E-step, and useMEX option
% 2017-08-17: ABC - added "full" covariance matrices
% 2017-11-25: ABC - small bug fix: invalid index when checking covar type
% 2018-02-13: v0.72 - check degenerate HMM emission has zero count
% 2018-05-24: v0.73 - minor bug fix: ZOmega
% 2018-04-24: v0.73 - add mu2 computation for MANOVA test (still testing)
% 2020-09-07: v0.78 - fix bug: ignore Z for entries with base omega=0


% v0.73 - check if stats are needed
COMPUTE_MU2 = mopt.computeStats;

tic;


h3m_r.LogLs = [];

num_iter = 0;

% number of components in the base and reduced mixtures
Kb = h3m_b.K;
Kr = h3m_r.K;

% number of states 
N = size(h3m_b.hmm{1}.A,1);
Nr = size(h3m_r.hmm{1}.A,1);

% number of mixture components in each state emission probability
M = h3m_b.hmm{1}.emit{1}.ncentres;

% length of the virtual sample sequences
T = mopt.tau;

% dimension of the emission variable
dim = h3m_b.hmm{1}.emit{1}.nin;

% number of virtual samples 
virtualSamples = mopt.Nv;

% correct
N_i = virtualSamples * h3m_b.omega * Kb;
N_i = N_i(:);

% regularization value
if isfield(mopt,'reg_cov')
    reg_cov = mopt.reg_cov;
else
    reg_cov = 0;
end

% minimum number of iterations
if ~isfield(mopt,'min_iter')
    mopt.min_iter = 0;
end

if ~isfield(mopt,'story')
    mopt.story = 0;
end

% timing
here_time = tic;
if isfield(mopt,'start_time')
    a_time = mopt.start_time;
else
    a_time = here_time;
end

% save the H3Ms in each iteration
if mopt.story
    story = {};
    time = [];
end
if mopt.story
    story{end+1} = h3m_r;
    time(end+1) = toc(a_time);
end

% regularize
for j = 1 : Kr
  for n = 1 : Nr  % bug fix (was N before)
    switch(h3m_r.hmm{j}.emit{n}.covar_type)
      case 'diag'        
        h3m_r.hmm{j}.emit{n}.covars = h3m_r.hmm{j}.emit{n}.covars + reg_cov;
      case 'full'
        h3m_r.hmm{j}.emit{n}.covars = h3m_r.hmm{j}.emit{n}.covars + reg_cov*eye(dim);
    end
  end
end

switch mopt.inf_norm
  case ''
    inf_norm = 1;
  case 'n'
    inf_norm = virtualSamples / Kb;
  case {'tn' 'nt'}
    inf_norm = T * virtualSamples / Kb;
  case {'t'}
    inf_norm = T;  
end

smooth =  mopt.smooth;

% 2017-08-08 - ABC - handle MEX 
if ~isfield(mopt, 'useMEX')
  mopt.useMEX = 1;
end

canuseMEX = 0;
if (mopt.useMEX)
  if (exist('hem_hmm_bwd_fwd_mex')==3)
    canuseMEX = 1;
    for n=1:length(h3m_b.hmm)
      for j=1:length(h3m_b.hmm{n}.emit)
        if h3m_b.hmm{n}.emit{j}.ncentres ~= 1
          canuseMEX = 0;
          warning('hem_hmm_bwd_fwd_mex: ncentres should be 1');
        end
        %if ~strcmp(h3m_b.hmm{n}.emit{j}.covar_type, 'diag')
        %  canuseMEX = 0;
        %  warning('hem_hmm_bwd_fwd_mex: only diag is supported');
        %end
        if canuseMEX == 0
          break;
        end
      end
      if canuseMEX == 0
        break
      end
    end
    
  else
    canuseMEX = 0;
    warning('hem_hmm_bwd_fwd_mex: MEX implementation not found');
  end
end

if (canuseMEX)
  % get maxN and maxN2
  maxN  = max(cellfun(@(x) length(x.prior), h3m_b.hmm));
  maxN2 = max(cellfun(@(x) length(x.prior), h3m_r.hmm));
  
else
  
  % make cache for faster E-step
  % v0.73 - put into its own function
  h3m_b_hmm_emit_cache = make_h3m_b_cache(h3m_b);
  
  % initialize arrays for E-step
  L_elbo           = zeros(Kb,Kr);
  nu_1             = cell(Kb,Kr);
  update_emit_pr   = cell(Kb,Kr);
  update_emit_mu   = cell(Kb,Kr);
  update_emit_M    = cell(Kb,Kr);
  sum_xi           = cell(Kb,Kr);
end


% start looping variational E step and M step
while 1

  %%%%%%%%%%%%%%%%%%%%
  %%%    E-step    %%%
  %%%%%%%%%%%%%%%%%%%%
    
  if (canuseMEX)
    %% MEX implementation
    
    % 2017-11-25: ABC: BUG FIX: n->1, j->1.  just need to check covariance type
    switch(h3m_b.hmm{1}.emit{1}.covar_type)
      case 'diag'
        [L_elbo, nu_1, update_emit_pr, update_emit_mu, update_emit_M, sum_xi] = ...
          hem_hmm_bwd_fwd_mex(h3m_b.hmm,h3m_r.hmm,T,smooth, maxN, maxN2);

      case 'full'
        % cache logdets and inverse covariances
        logdetCovR = cell(1,Kr);
        invCovR    = cell(1,Kr);
        for j=1:Kr
           % does not support ncentres>1
          logdetCovR{j} = zeros(1,length(h3m_r.hmm{j}.prior));
          invCovR{j}    = zeros(dim, dim, length(h3m_r.hmm{j}.prior));
          for n=1:length(h3m_r.hmm{j}.prior)
            logdetCovR{j}(n)  = log(det(h3m_r.hmm{j}.emit{n}.covars)); 
            invCovR{j}(:,:,n) = inv(h3m_r.hmm{j}.emit{n}.covars);
          end
        end
        
        [L_elbo, nu_1, update_emit_pr, update_emit_mu, update_emit_M, sum_xi] = ...
          hem_hmm_bwd_fwd_mex(h3m_b.hmm,h3m_r.hmm,T,smooth, maxN, maxN2, logdetCovR, invCovR);
      
      otherwise  
        error('not supported')
    end
    
    % consistency check
    if 0
      err = [];
      err(end+1) = totalerror(xL_elbo, L_elbo);
      err(end+1) = totalerror(xnu_1, nu_1);
      err(end+1) = totalerror(xsum_xi, sum_xi);
      err(end+1) = totalerror(xupdate_emit_pr, update_emit_pr);
      err(end+1) = totalerror(xupdate_emit_mu, update_emit_mu);
      err(end+1) = totalerror(xupdate_emit_M, update_emit_M);
      if any(err>1e-6)
        warning('mismatch');
        keyboard
      end
    end
    
  else
    %% MATLAB implementation
    % loop reduced mixture components
    for j = 1 : Kr
      
      % get this HMM in the reduced mixture
      hmm_r = h3m_r.hmm{j};
      
      % loop base mixture components
      for i = 1 : Kb
        
        % get this HMM in the base mixture
        hmm_b = h3m_b.hmm{i};
        
        % run the Hierarchical bwd and fwd recursions. This will compute:
        % the ELBO on E_M(b)i [log P(y_1:T | M(r)i)]
        % the nu_1^{i,j} (rho)
        % the sum_gamma  [sum_t nu_t^{i,j} (rho,gamma)] sum_m c.. nu..
        % the sum_gamma  [sum_t nu_t^{i,j} (rho,gamma)] sum_m c.. nu.. mu ..
        % the sum_gamma  [sum_t nu_t^{i,j} (rho,gamma)] sum_m c.. nu.. M .. which accounts for the mu mu' and the sigma
        % [the last 3 are not normalized, of course]
        % the sum_t xi_t^{i,j} (rho,sigma)     (t = 2 ... T)
        [L_elbo(i,j) ...
          nu_1{i,j} ...
          update_emit_pr{i,j} ...
          update_emit_mu{i,j} ...
          update_emit_M{i,j}  ...
          sum_xi{i,j} ...
          ] = hem_hmm_bwd_fwd(hmm_b,hmm_r,T,smooth, h3m_b_hmm_emit_cache{i});
        
        % consistency check
        if 0
          [xL_elbo ...
            xnu_1 ...
            xupdate_emit_pr ...
            xupdate_emit_mu ...
            xupdate_emit_M  ...
            xsum_xi ...
            ] = hem_hmm_bwd_fwd_OLD(hmm_b,hmm_r,T,smooth);
          err = [];
          err(end+1) = abs(xL_elbo-L_elbo(i,j));
          err(end+1) = sum(abs(xnu_1(:) - nu_1{i,j}(:)));
          err(end+1) = sum(abs(xupdate_emit_pr(:) - update_emit_pr{i,j}(:)));
          err(end+1) = sum(abs(xupdate_emit_mu(:) - update_emit_mu{i,j}(:)));
          err(end+1) = sum(abs(xupdate_emit_M(:)  - update_emit_M{i,j}(:)));
          err(end+1) = sum(abs(xsum_xi(:) - sum_xi{i,j}(:)));
          if any(err>1e-6)
            warning('mismatch')
            keyboard;
          end
          
        end
        % end loop base
      end
      
      % end loop reduced
    end
  end
  
    
  % compute the z_ij
  % this is not normalized ...
  % log_Z = ones(Kb,1) * log(h3m_r.omega) + diag(N_i) * L_elbo;
  L_elbo = L_elbo / (inf_norm);
  log_Z = ones(Kb,1) * log(h3m_r.omega) + (N_i(:) * ones(1,Kr)) .* L_elbo;
    
  % normalize Z
  Z = exp(log_Z - logtrick(log_Z')' * ones(1,Kr));
  
  % compute the elbo to total likelihood ...  
  %     ll = Z .* (log_Z - log(Z));
  %
  %     ll(isnan(ll)) = 0;
  %
  %     % new_LogLikelihood = sum(sum(ll));
  %     new_LogLikelihood_foo = ones(1,Kb) * ll * ones(Kr,1);
    
  % these should be the same
  new_LogLikelihood = sum(logtrick(log_Z')');
    
  % update the log likelihood in the reduced mixture
  old_LogLikelihood = h3m_r.LogL;
  h3m_r.LogL = new_LogLikelihood;
  h3m_r.LogLs(end+1) = new_LogLikelihood;
  h3m_r.Z = Z;
    
  % check whether to continue with a new iteration or to return
  stop = 0;
  if num_iter > 1
    changeLL = (new_LogLikelihood - old_LogLikelihood) / abs(old_LogLikelihood);
  else
    changeLL = inf;
  end
  
  if (mopt.verbose > 2)
    fprintf('iter=%d; LL=%g (pL=%g)\n', num_iter, new_LogLikelihood, changeLL);
  end
  
  if changeLL<0
    warning('The change in log likelihood is negative!!!')
  end
  
  % change in LL is below threshold?
  switch mopt.termmode
    case 'L'
      if (changeLL < mopt.termvalue)
        stop = 1;
      end
  end
  
  % max iteration is reached?
  if (num_iter > mopt.max_iter)
    stop = 1;
  end
        
  % stop and reached minimum number of iterations?
  % then exit the EM loop
  if (stop) && (num_iter >= mopt.min_iter)
    
    if COMPUTE_MU2
      %% v0.73 - converged, now compute mu2 for MANOVA test
      
      % TODO: this is not very efficient because it recomputes everything.
      % we can cache the values from E-step, and compute just the mu2 term here.
      
      % Caching/reshaping for faster E-step
      h3m_b_hmm_emit_cache = make_h3m_b_cache(h3m_b);
      xL_elbo           = zeros(Kb,Kr);
      xnu_1             = cell(Kb,Kr);
      xupdate_emit_pr   = cell(Kb,Kr);
      xupdate_emit_mu   = cell(Kb,Kr);
      xupdate_emit_M    = cell(Kb,Kr);
      xsum_xi           = cell(Kb,Kr);
      xupdate_emit_mu2  = cell(Kb,Kr);
      
      % Run E-step again to get mu2
      % loop reduced mixture components
      for j = 1 : Kr
        % get this HMM in the reduced mixture
        hmm_r = h3m_r.hmm{j};
        % loop base mixture components
        for i = 1 : Kb
          % get this HMM in the base mixture
          hmm_b = h3m_b.hmm{i};
          
          % get mu2 - USE MEX VERSION (slow but only need to do it once)
          [xL_elbo(i,j), ...
            xnu_1{i,j}, ...
            xupdate_emit_pr{i,j}, ...
            xupdate_emit_mu{i,j}, ...
            xupdate_emit_M{i,j},  ...
            xsum_xi{i,j}, ...
            ~, ...
            xupdate_emit_mu2{i,j} ...
            ] = hem_hmm_bwd_fwd(hmm_b,hmm_r,T,smooth, h3m_b_hmm_emit_cache{i});
        end
      end
      
      % for each reduced component j
      for j = 1:Kr
        xestats.Zomega   = Z(:,j) .* h3m_b.omega(:);
        xestats.nu_1     = xnu_1(:,j);
        xestats.sum_xi   = xsum_xi(:,j);
        xestats.emit_pr  = xupdate_emit_pr(:,j);
        xestats.emit_mu  = xupdate_emit_mu(:,j);
        xestats.emit_M   = xupdate_emit_M(:,j);
        xestats.emit_mu2 = xupdate_emit_mu2(:,j);
        
        % setup constants
        xestats.M   = M;
        xestats.N2  = size(h3m_r.hmm{j}.A,1);
        xestats.dim = dim;
        xestats.Kb  = Kb;
        xestats.reg_cov = reg_cov;
        
        % run M-step on the j-th component
        tmphmm = hem_mstep_component(h3m_r.hmm{j}, xestats, mopt);
        
        % copy statistics
        h3m_r.hmm{j}.stats = tmphmm.stats;
      end
    end
    
    %% exit the EM loop
    break
  end
  
  num_iter = num_iter + 1;  
  old_LogLikelihood = new_LogLikelihood;
    
  %%%%%%%%%%%%%%%%%%%%
  %%%    M-step    %%%
  %%%%%%%%%%%%%%%%%%%%
  % compute new parameters

  % make a copy first
  h3m_r_new = h3m_r;
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%  re-estimation of the component weights (omega)
  
  % v0.78 - ABC - 2020-09-07
  % BUG FIX: ignore Z values where h3m_b.omega has 0 entries
  %     this only happens when doing initialization for co-clustering,
  %     and there are missing HMMs for some stimuli
  
  %omega_new = (ones(1,Kb) / Kb) * Z;
  ignoreZ = sign(h3m_b.omega(:));
  omega_new = sum(Z.*ignoreZ, 1) / sum(ignoreZ);
  h3m_r_new.omega = omega_new;
  
  %% for each reduced component j
  for j = 1:Kr    
    % get statistics for j-th component
    %estats.Zomega  = Z(:,j) .* N_i(:);     
    % BUG FIX: v0.73 - N_i contains an extra factor of Nvs*Kb.
    % it doesn't affect the M-step, since the extra factors cancel each other.
    % It does affect count_emit, which has an extra multiplication factor.
    % Previous usage of count_emit, only foudn the max index or equal-to-zero,
    % so the bug has no effect.
    % However, now we want to use count_emit for statistical tests, so we need the actual number.
    estats.Zomega  = Z(:,j) .* h3m_b.omega(:);
    estats.nu_1    = nu_1(:,j);
    estats.sum_xi  = sum_xi(:,j);
    estats.emit_pr = update_emit_pr(:,j);
    estats.emit_mu = update_emit_mu(:,j);
    estats.emit_M  = update_emit_M(:,j);
    
    % setup constants
    estats.M   = M;
    estats.N2  = size(h3m_r.hmm{j}.A,1);
    estats.dim = dim;
    estats.Kb  = Kb;
    estats.reg_cov = reg_cov;
    
    % run M-step on the j-th component
    h3m_r_new.hmm{j} = hem_mstep_component(h3m_r_new.hmm{j}, estats, mopt);
  end
    
  %% check for degenerate components %%%%%%%%%%%%%%%%%
        
  % if an HMM component gets 0 prior ...
  ind_zero = find(h3m_r_new.omega == 0);
  for i_z = ind_zero;
    h3m_r_new = hem_fix_degenerate_component(h3m_r_new, i_z, N, Nr);
  end
  
  % ABC: 2018-02-13
  % check if an HMM state has zero counts
  for j = 1 : Kr
    ind_zero = find(h3m_r_new.hmm{j}.stats.emit_vcounts == 0);
    if ~isempty(ind_zero)
      for i_z = ind_zero
        h3m_r_new.hmm{j} = hem_fix_degenerate_hmm(h3m_r_new.hmm{j}, i_z);
      end
    end
  end
  
  
  % check if an emission Gaussian has zero prior (only applicable for GMM emission with >1 components
  if (h3m_r_new.hmm{1}.emit{1}.ncentres > 1)
    for j = 1 : Kr
      N2 = size(h3m_r.hmm{j}.A,1);      
      for n = 1 : N2
        % if some of the emission prob is zero, replace it ...
        ind_zero = find(h3m_r_new.hmm{j}.emit{n}.priors == 0);
        for i_z = ind_zero;
          h3m_r_new.hmm{j}.emit{n} = hem_fix_degenerate_emission(h3m_r_new.hmm{j}.emit{n}, i_z);
        end
      end
    end
  end
  
        
  %% update the model %%%%%%%%%%%%%%%%%%%%
  h3m_r = h3m_r_new;
  
  
  % add to story
  if mopt.story
    story{end+1} = h3m_r;
    time(end+1) = toc(a_time);
  end  
end


h3m_r.elapsed_time = toc(here_time);
if mopt.story
  h3m_r.story = story;
  h3m_r.time = time;
end

% end of the function
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s] = logtrick(lA)
% logtrick - "log sum trick" - calculate log(sum(A)) using only log(A) 
%
%   s = logtrick(lA)
%
%   lA = column vector of log values
%
%   if lA is a matrix, then the log sum is calculated over each column
% 

[mv, mi] = max(lA, [], 1);
temp = lA - repmat(mv, size(lA,1), 1);
cterm = sum(exp(temp),1);
s = mv + log(cterm);
end



%% make cache for non-MEX version of E-step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h3m_b_hmm_emit_cache = make_h3m_b_cache(h3m_b)

Kb = h3m_b.K;

% 2017-08-07 - ABC - Caching/reshaping for faster E-step
h3m_b_hmm_emit_cache = cell(1,Kb);
for i=1:Kb
  g3m_b = h3m_b.hmm{i}.emit;
  BKall = cellfun(@(x) x.ncentres, g3m_b);
  BKmax = max(BKall);
  
  % does not support different number of centres (for speed)
  if ~all(BKall==BKmax)
    error('different number of centres is not supported.');
  end
  
  gmmBcentres_tmp = cellfun(@(x) x.centres, g3m_b, ...
    'UniformOutput', false);
  h3m_b_hmm_emit_cache{i}.gmmBcentres = cat(3,gmmBcentres_tmp{:});  % BKmax x dim x N
  
  % extract all the covars
  gmmBcovars_tmp = cellfun(@(x) x.covars, g3m_b, ...
    'UniformOutput', false);
  switch (g3m_b{1}.covar_type)
    case 'diag'
      h3m_b_hmm_emit_cache{i}.gmmBcovars = cat(3,gmmBcovars_tmp{:});  % BKmax x dim x N
    case 'full'
      h3m_b_hmm_emit_cache{i}.gmmBcovars = cat(4,gmmBcovars_tmp{:});  % dim x dim x BKmax x N
      % TODO: add cache for logdet(gmmBcovars)
  end
  
  % extract all priors
  gmmBpriors_tmp = cellfun(@(x) x.priors', g3m_b, ...
    'UniformOutput', false);
  h3m_b_hmm_emit_cache{i}.gmmBpriors = cat(2,gmmBpriors_tmp{:});  % BKmax x N
end

end



