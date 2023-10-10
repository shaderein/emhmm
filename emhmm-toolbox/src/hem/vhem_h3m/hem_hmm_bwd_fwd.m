function [LL_elbo, ...
          sum_nu_1, ...
          update_emit_pr, ...
          update_emit_mu, ...
          update_emit_Mu,  ...
          sum_xi, ...
          sum_t_nu, ...
          update_emit_mu2] = hem_hmm_bwd_fwd(hmm_b,hmm_r,T,smooth, h3m_b_hmm_emit_cache)
% run the Hierarchical bwd and fwd recursions. This will compute:
% the ELBO on E_M(b)i [log P(y_1:T | M(r)i)]
% the nu_1^{i,j} (rho)
% the sum_gamma  [sum_t nu_t^{i,j} (rho,gamma)] sum_m c.. nu..
% the sum_gamma  [sum_t nu_t^{i,j} (rho,gamma)] sum_m c.. nu.. mu ..
% the sum_gamma  [sum_t nu_t^{i,j} (rho,gamma)] sum_m c.. nu.. M .. which accounts for the mu mu' and the sigma
% [the last 3 are NOT normalized YET, of course]
% the sum_t xi_t^{i,j} (rho,sigma)     (t = 2 ... T)
% ---
% H3M Toolbox 
% Copyright 2012, Emanuele Coviello, CALab, UCSD

% 2017-08-06 - ABC - MATLAB speed optimziations 
% 2017-08-18 - ABC - added full covariance matrices
% 2018-05-25 - v0.73 - added update_emit_mu2 for MANOVA test

if isempty(smooth)
    smooth = 1;
end

% v0.73 - compute mu2?
if nargout < 8
  COMPUTE_MU2 = 0;
else
  COMPUTE_MU2 = 1;
end

% number of states 
N = size(hmm_b.A,1);
N2 = size(hmm_r.A,1);
% number of mixture components in each state emission probability
M = hmm_b.emit{1}.ncentres;
% dimension of the emission variable
dim = hmm_b.emit{1}.nin;

% first, get the elbo of the E_gauss(b) [log p (y | gauss (r))]  (different b on different rows)
% and the sufficient statistics for the later updates
% sum_w_pr is a N by 1 cell of N by M
% sum_w_mu is a N by 1 cell of N by dim by M
% sum_w_Mu is a N by 1 cell of N by dim by M (for the diagonal case)
% sum_w_Mu is a N by 1 cell of N by dim x dim by M (for the full case)
% sum_w_mu2 is a N by 1 cell of N by dim by dim by M

if COMPUTE_MU2
  % v0.73 - also compute mu2
  [LLG_elbo, sum_w_pr, sum_w_mu, sum_w_Mu, sum_w_mu2] = g3m_stats(hmm_b.emit,hmm_r.emit, h3m_b_hmm_emit_cache);
else
  [LLG_elbo, sum_w_pr, sum_w_mu, sum_w_Mu] = g3m_stats(hmm_b.emit,hmm_r.emit, h3m_b_hmm_emit_cache);
end


% consistency check
if 0
  [LLG_elbo_OLD sum_w_pr_OLD sum_w_mu_OLD sum_w_Mu_OLD] = g3m_stats_OLD(hmm_b.emit,hmm_r.emit);
  
  errs = [];
  errs(end+1) = sum(sum(abs(LLG_elbo_OLD-LLG_elbo)));
  errs(end+1) = sum(sum(sum(sum(abs(cat(4, sum_w_pr_OLD{:})-cat(4, sum_w_pr{:}))))));
  errs(end+1) = sum(sum(sum(sum(abs(cat(4, sum_w_mu_OLD{:})-cat(4, sum_w_mu{:}))))));
  errs(end+1) = sum(sum(sum(sum(abs(cat(4, sum_w_Mu_OLD{:})-cat(4, sum_w_Mu{:}))))));
  
  if any(errs>1e-10)
    warning('mismatch');
    keyboard
  end
end

LLG_elbo = LLG_elbo / smooth; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% do the backward recursion %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ab = hmm_b.A;
Ar = hmm_r.A;

% allocate Theta, i.e., the VA parameters in the form of a CMC
Theta = zeros(N2,N2,N,T); % the dimensions refer to rho sigma (states of M(r)) gamma (state of M(b)) and time
        % rho sigma gamma t

% allocate the log-likelihood
LL_old = zeros(N,N2);       % the dimensions refer to gamma (state of M(b)) and sigma (states of M(r)) 

% the log-likelihood can be intialized by all zeros ...

for t = T : -1 : 2
    
  % [N,N2]
  LL_new = zeros(size(LL_old));

  if 0
    %% old code
    for rho = 1 : N2
      
      % [N2xN] =          [1xN2]' x [1xN]      + [N,N2]'   +  [N,N2]'
      logtheta = (log(Ar(rho,:))' * ones(1,N)) + LLG_elbo' + LL_old';
      % % % logtheta = (log(Ar(rho,:))' * ones(1,N)) + LLG_elbo + LL_old;
      
      % [1xN]
      logsumtheta = logtrick(logtheta);
      
      % [Nx1]       =  [NxN] * [1xN]'
      LL_new(:,rho) =  Ab * logsumtheta';
      
      % normalize so that each clmn sums to 1 (may be not necessary ...)
      % [N2xN] = [N2xN]   -       [N2x1]*[1xN]
      theta = exp(logtheta - ones(N2,1) * logsumtheta);
      
      % and store for later
      % [1xN2xNx1]
      Theta(rho,:,:,t) = theta;
      
    end
  end
  
  %% 2017-08-07 - ABC - fast code
  % [N2xN]xN2 =          [1xN2]' x [1xN]      + [N,N2]'   +  [N,N2]'
  %    logtheta = (log(Ar(rho,:))' * ones(1,N)) + LLG_elbo' + LL_old';
  logtheta_all = bsxfun(@plus, ...
    reshape(log(Ar)', [N2 1 N2]), ...
    reshape(LLG_elbo',[N2 N]) + reshape(LL_old',[N2 N]));
 
  % [1xN]xN2
  %logsumtheta = logtrick(logtheta);
  logsumtheta_all = logtrick2(logtheta_all, 1);
    
  % [NxN2]         =  [NxN] * [1xN]'
  %  LL_new(:,rho) =  Ab * logsumtheta';
  LL_new = reshape(sum(bsxfun(@times, reshape(Ab,[N N 1]), logsumtheta_all), 2), [N N2]);
  
  % normalize so that each clmn sums to 1 (may be not necessary ...)
  % [N2xN]xN2 = [N2xN]   -       [N2x1]*[1xN]
  %theta = exp(logtheta - ones(N2,1) * logsumtheta);
  theta_all = exp(bsxfun(@minus, logtheta_all, logsumtheta_all));
  
  % and store for later
  % [1xN2xNx1]
  %Theta(rho,:,:,t) = theta;
  Theta(:,:,:,t) = permute(theta_all, [3 1 2]);
  
  %% consistency check
  if 0
    err = [];
    err(end+1) = sum(sum(abs(LL_new - xLL_new)));
    err(end+1) = sum(sum(sum(abs(Theta(:,:,:,t) - xTheta(:,:,:,t)))));
    if any(err>1e-10)
      warning('mismatch error');
      keyboard
    end
  end
  
  LL_old = LL_new;
end

% terminate the recursion

logtheta = (log(hmm_r.prior) * ones(1,N)) + LLG_elbo' + LL_old';
% % % logtheta = (log(hmm_r.prior) * ones(1,N)) + LLG_elbo + LL_old;

logsumtheta = logtrick(logtheta);

LL_elbo =  hmm_b.prior' * logsumtheta';

% normalize so that each clmn sums to 1 (may be not necessary ...) 
theta = exp(logtheta - ones(N2,1) * logsumtheta);

% and store for later
Theta_1 = theta;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% do the forward recursion  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rather then saving all intermediate values, just cache the cumulative
% values that are needed for the updates (it saves a lot of memory)

% nu is N by N, first dimension indexed by sigma (M(r)), second by gamma (M(b))

% initialize [N2 x N]
nu = (ones(N2,1) * hmm_b.prior') .* Theta_1; % N by N (indexed by sigma and gamma)



% CACHE: sum_gamma nu_1(sigma,gamma) (this is one of the outputs ...)
sum_nu_1 = sum(nu,2)';

% CACHE: sum_t nu_t(sigma, gamma)
% sum_t_nu = zeros(N,N);
sum_t_nu = nu;

% CACHE: sum_t sum_gamma xi(rho,sigma,gamma,t)
sum_t_sum_g_xi = zeros(N2,N2); % N by N (indexed by rho and sigma)

for t = 2 : T
    
    % compute the inner part of the update of xi (does not depend on sigma)
    foo = nu * Ab; % indexed by rho gamma [N2xN]
    
    %% old code using loop
    if 0
      for sigma = 1 : N2
        % new xi
        % xi(:,sigma,:,t) = foo .* squeeze(Theta(:,sigma,:,t));
        %xi_foo = foo .* squeeze(Theta(:,sigma,:,t)); % (indexed by rho gamma);
        
        % ABC: bug fix when another dim is 1
        xi_foo = foo .* reshape(Theta(:,sigma,:,t), [size(Theta,1), size(Theta,3)]); % (indexed by rho gamma);
        
        % CACHE:
        sum_t_sum_g_xi(:,sigma) = sum_t_sum_g_xi(:,sigma) + sum(xi_foo,2);
        
        % new nu
        % nu(sigma,:) = ones(1,N) * squeeze(xi(:,sigma,:,t));
        nu(sigma,:) = ones(1,N2) * xi_foo;
      end
    end
    
    %% 2018-08-07 - ABC - faster code w/o loop
    % xi_foo = foo .* reshape(Theta(:,sigma,:,t), [size(Theta,1), size(Theta,3)]); % (indexed by rho gamma);
    xi_foo_all = bsxfun(@times, reshape(foo,[N2 1 N]), Theta(:,:,:,t));
    
    % sum_t_sum_g_xi(:,sigma) = sum_t_sum_g_xi(:,sigma) + sum(xi_foo,2);
    sum_t_sum_g_xi = sum_t_sum_g_xi + sum(xi_foo_all, 3);
    
    % nu(sigma,:) = ones(1,N2) * xi_foo;
    nu = reshape(sum(xi_foo_all,1), [N2 N]);
    
    %% consistency check
    if 0
      err = [];
      err(end+1) = sum(sum(abs(xnu-xnu)));
      err(end+1) = sum(sum(abs(xsum_t_sum_g_xi - sum_t_sum_g_xi)));
      if any(err>1e-10)
        warning('mismatch error');
        keyboard
      end
    end
    
    % CACHE: in the sum_t nu_t(sigma, gamma)
    sum_t_nu = sum_t_nu + nu;
    
end

% this is one of the outputs ...
sum_xi = sum_t_sum_g_xi;


%%%% now prepare the cumulative sufficient statistics for the reestimation
%%%% of the emission distributions

update_emit_pr = zeros(N2,M);
update_emit_mu = zeros(N2,dim,M);
switch hmm_b.emit{1}.covar_type
  case 'diag'
    update_emit_Mu = zeros(N2,dim,M);
  case 'full'
    update_emit_Mu = zeros(N2,dim,dim,M);
end

% v0.73 - initialize mu2
if COMPUTE_MU2
  update_emit_mu2 = zeros(N2,dim,dim,M);
end  

% loop all the emission GMM of each state
for sigma = 1 : N2
  
  update_emit_pr(sigma,:) = sum_t_nu(sigma,:) * sum_w_pr{sigma};
  
  foo_sum_w_mu = sum_w_mu{sigma};
  foo_sum_w_Mu = sum_w_Mu{sigma};
  
  for l = 1 : M
    update_emit_mu(sigma,:,l) = sum_t_nu(sigma,:) * foo_sum_w_mu(:,:,l);
    
    switch hmm_b.emit{1}.covar_type
      case 'diag'
        update_emit_Mu(sigma,:,l) = sum_t_nu(sigma,:) * foo_sum_w_Mu(:,:,l);
      case 'full'
        update_emit_Mu(sigma,:,:,l) = sum(bsxfun(@times,  ...
          reshape(sum_t_nu(sigma,:), [N 1 1]), foo_sum_w_Mu(:,:,:,l)), 1);
    end
  end
  
  % v0.73 - compute mu2
  if COMPUTE_MU2
    foo_sum_w_mu2 = sum_w_mu2{sigma};
    for k = 1:M
      update_emit_mu2(sigma,:,:,l) = sum(bsxfun(@times,  ...
        reshape(sum_t_nu(sigma,:), [N 1 1]), foo_sum_w_mu2(:,:,:,l)), 1);
    end 
  end
  
end


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s] = logtrick2(lA, dim)
% logtrick - "log sum trick" - calculate log(sum(A)) using only log(A) 
%
%   s = logtrick(lA, dim)
%
%   lA = matrix of log values
%
%   if lA is a matrix, then the log sum is calculated over dimension dim.
% 

[mv, mi] = max(lA, [], dim);
rr = size(lA);
rr(dim) = 1;
temp = bsxfun(@minus, lA, reshape(mv, rr));
cterm = sum(exp(temp),dim);
s = mv + log(cterm);
end



 