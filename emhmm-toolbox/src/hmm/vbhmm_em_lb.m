function [LB, d_LB] = vbhmm_em_lb(lbhmm_stats, fbhmm_varpar, fbstats)
% vbhmm_em_lb - calculate log-likelihood lower-bound
%
% this is an internal function called by vbhmm_em
%
% [LB] = vbhmm_em_lb(lbhmm_stats, hmm_varpar, fbstats)
%
% [LB, d_LB] = vbhmm_em_lb(lbhmm_stats, hmm_varpar, fbstats)
%  set lbhmm_stats.do_deriv = 1;
%  calculate the derivatives wrt the hyperparameters
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-01-13
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong

% 2017-06-09: ABC - added option to calculate derivatives of LB
% 2017-06-14: ABC - epsilon & logATilde changed to row format
% 2017-07-31: ABC - clean up code for derivatives
% 2017-08-21: ABC - set derivatives to zero when hyps are clipped at extreme values

% 2018-02-11: v0.72 - xxx remove hyp clipping
%             v0.72 - compute log(det(W0inv)) more accurately for large W0
% 2018-03-10: v0.72 - fixed bugs with group_inds not passed
% 2018-05-15: v0.72 - use a loop over clipped.

%% extract the constants, stats, etc

% hyperparameters
dim       = lbhmm_stats.dim;
K         = lbhmm_stats.K;
N         = lbhmm_stats.N;
alpha0    = lbhmm_stats.alpha0;
epsilon0  = lbhmm_stats.epsilon0;
m0        = lbhmm_stats.m0;
beta0     = lbhmm_stats.beta0;
v0        = lbhmm_stats.v0;
W0inv     = lbhmm_stats.W0inv;
W0mode    = lbhmm_stats.W0mode;
usegroups = lbhmm_stats.usegroups;
if usegroups
  group_inds = lbhmm_stats.group_inds;
end
t1_S      = lbhmm_stats.t1_S;
xbar      = lbhmm_stats.xbar;
Nk        = lbhmm_stats.Nk;
M         = lbhmm_stats.M;
clipped   = lbhmm_stats.clipped;

if isfield(lbhmm_stats, 'do_deriv')
  do_deriv = lbhmm_stats.do_deriv;
else
  do_deriv = 0;
end

% hmm parameters: v, W, epsilon, alpha, m, beta
v       = fbhmm_varpar.v;
W       = fbhmm_varpar.W;
epsilon = fbhmm_varpar.epsilon;
alpha   = fbhmm_varpar.alpha;
m       = fbhmm_varpar.m;
beta    = fbhmm_varpar.beta;

% get statistics and constants
logrho_Saved  = fbstats.logrho_Saved;
t_gamma_Saved = fbstats.gamma_all;
xi_sum = fbstats.xi_sum;
phi_norm  = fbstats.phi_norm;
logLambdaTilde = fbstats.logLambdaTilde;
logPiTilde = fbstats.logPiTilde;
logATilde = fbstats.logATilde;

%% calculate constants
switch(W0mode)
  case 'iid'
    logdetW0inv = dim*log(W0inv(1,1));    
  case 'diag'
    logdetW0inv = sum(log(diag(W0inv)));
  otherwise
    error('bad W0mode');
end

logCalpha0 = gammaln(K*alpha0) - K*gammaln(alpha0);
for k = 1:K
  logCepsilon0(k) = gammaln(K*epsilon0) - K*gammaln(epsilon0);
end
logB0 = (v0/2)*logdetW0inv - (v0*dim/2)*log(2) ...
  - (dim*(dim-1)/4)*log(pi) - sum(gammaln(0.5*(v0+1-[1:dim])));

if ~usegroups
  logCalpha = gammaln(sum(alpha)) - sum(gammaln(alpha));
  for k = 1:K
    logCepsilon(k) = gammaln(sum(epsilon(k,:))) - sum(gammaln(epsilon(k,:)));
  end
else
  numgroups = length(alpha);
  for g=1:numgroups
    logCalpha{g} = gammaln(sum(alpha{g})) - sum(gammaln(alpha{g}));
    for k = 1:K
      logCepsilon{g}(k) = gammaln(sum(epsilon{g}(k,:))) - sum(gammaln(epsilon{g}(k,:)));
    end
  end
end
    
H = 0;
for k = 1:K
  logBk = -(v(k)/2)*log(det(W(:,:,k))) - (v(k)*dim/2)*log(2)...
    - (dim*(dim-1)/4)*log(pi) - sum(gammaln(0.5*(v(k) + 1 - [1:dim])));
  H = H - logBk - 0.5*(v(k) - dim - 1)*logLambdaTilde(k) + 0.5*v(k)*dim;
  trSW(k) = trace(t1_S(:,:,k)*W(:,:,k));
  xbarT = xbar(k,:)';
  diff = xbarT - m(:,k);
  xbarWxbar(k) = diff'*W(:,:,k)*diff;
  diff = m(:,k) - m0;
  mWm(k) = diff'*W(:,:,k)*diff;
  trW0invW(k) = trace(W0inv*W(:,:,k));
end
    
%% calculate each term in the LB

% E(log p(X|Z,mu,Lambda)  Bishop (10.71) - ABC term 1
Lt1 = 0.5*sum(Nk.*(logLambdaTilde' - dim./beta...
  - v.*trSW' - v.*xbarWxbar' - dim*log(2*pi)));
    
% initial responsibilities (t=1)
gamma1 = t_gamma_Saved(:,:,1);
    
% transition responsibilities
%gamma_t1 = zeros(K,K,N);
%for i = 1:size(logPiTilde,1)
%  for k = 1:size(gamma1,1)
%    for n = 1:size(gamma1,2)
%      gamma_t1(i,k,n) = gamma1(k,n).*logPiTilde(i);
%    end
%  end
%end
    
% E[log p(Z|pi)]   Bishop (10.72) - ABC term 2, part 1
if ~usegroups
  PiTilde_t = repmat(logPiTilde,1,N);
  gamma_t1 = gamma1.*PiTilde_t;
  Lt2a = sum(sum(gamma_t1));
else
  Lt2a = 0;
  for g=1:numgroups
    PiTilde_t = logPiTilde{g};
    gamma_t1 = bsxfun(@times, gamma1(:,group_inds{g}), PiTilde_t);
    Lt2a = Lt2a + sum(sum(gamma_t1));
  end
end
      
% E[log p(Z|A)]   ~Bishop 10.72 - ABC term 2, part 2 [CORRECT?]
if ~usegroups
  %sumxi_t1 = zeros(N,K,K);
  %for i = 1:size(t_sumxi_Saved,1)
  %  sumxi_t1(i,:,:) = squeeze(t_sumxi_Saved(i,:,:)).*ATilde_t;
  %end
  %Lt2b = sum(sum(sum(sumxi_t1)))
  Lt2b = sum(M(:).*logATilde(:));
else
  Lt2b = 0;
  for g=1:numgroups
    Lt2b = Lt2b + sum(M{g}(:).*logATilde{g}(:));
  end
end

% E[log p(Z|pi, A)]
% ABC term 2
Lt2 = Lt2a + Lt2b;

% E[log p(pi)]   Bishop (10.73)   ABC term 3
if ~usegroups
  Lt3 = logCalpha0 + (alpha0-1)*sum(logPiTilde);
else
  Lt3 = 0;
  for g=1:numgroups
    Lt3 = Lt3 + logCalpha0 + (alpha0-1)*sum(logPiTilde{g});
  end
end

% E[log p(A)] = sum E[log p(a_j)]   (equivalent to Bishop 10.73) ABC term 4
if ~usegroups
  for k = 1:K
    Lt4a(k) = logCepsilon0(k) + (epsilon0 -1)*sum(logATilde(k,:));
  end
  Lt4 = sum(Lt4a);
else
  Lt4 = 0;
  for g=1:numgroups
    for k = 1:K
      Lt4a(k) = logCepsilon0(k) + (epsilon0 -1)*sum(logATilde{g}(k,:));
    end
    Lt4 = Lt4 + sum(Lt4a);
  end
end

% E[log p(mu, Lambda)]  Bishop (10.74)  ABC term 5
Lt51 = 0.5*sum(dim*log(beta0/(2*pi)) + logLambdaTilde' - dim*beta0./beta - beta0.*v.*mWm');
Lt52 = K*logB0 + 0.5*(v0-dim-1)*sum(logLambdaTilde) - 0.5*sum(v.*trW0invW');
Lt5 = Lt51+Lt52;

% 2016-04-26 ABC: use correct q(Z)
    
% 2016-04-29 ABC: E[z log pi] (same as Lt2a)
%Lt61 = sum(sum(bsxfun(@times, t_gamma_Saved(:,:,1), logPiTilde)));
Lt61 = Lt2a;
    
% 2016-04-29 ABC:  E[zt zt-1 log a]  (same as Lt2b)
Lt62 = Lt2b;

% 2016-04-29 ABC:  E[z log rho]
% the zeros in logrho should remove times (tT+1):N
Lt63 = sum(sum(sum(t_gamma_Saved.*logrho_Saved)));

% 2016-04-29 ABC: normalization constant for q(Z)
Lt64 = sum(phi_norm);
%fprintf('   norm constant: %g\n', Lt64);
    
% 2016-04-29 ABC: E[log q(Z)] - ABC term 6
Lt6 = Lt61 + Lt62 + Lt63 - Lt64;
    
% E[log q(pi)]  Bishop (10.76)
if ~usegroups
  Lt71 = sum((alpha - 1).*logPiTilde) + logCalpha;
else
  Lt71 = 0;
  for g=1:numgroups
    Lt71 = Lt71 + sum((alpha{g} - 1).*logPiTilde{g}) + logCalpha{g};
  end
end

% E[log q(aj)]  (equivalent to Bishop 10.76)
if ~usegroups
  for k = 1:K
    Lt72(k) = sum((epsilon(k,:) - 1).*logATilde(k,:)) + logCepsilon(k);
  end
  Lt72sum = sum(Lt72);
else
  Lt72sum = 0;
  for g=1:numgroups
    for k = 1:K
      Lt72(k) = sum((epsilon{g}(k,:) - 1).*logATilde{g}(k,:)) + logCepsilon{g}(k);
    end
    Lt72sum = Lt72sum + sum(Lt72);
  end
end

% E[log q(pi, A)] - ABC term 7
Lt7 = Lt71 + Lt72sum;

% E[q(mu,Lamba)]  Bishop (10.77) - ABC term 8
Lt8 = 0.5*sum(logLambdaTilde' + dim.*log(beta/(2*pi))) - 0.5*dim*K - H;

%% Lower-bound value
% sum all terms together
LB = Lt1 + Lt2 + Lt3 + Lt4 + Lt5 - Lt6 - Lt7 - Lt8;


%% Calculate the derivatives %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (do_deriv)
  
  %% alpha0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  d_logC0.d_alpha0 = K*psi(0,K*alpha0) - K*psi(0,alpha0);
  if ~usegroups
    sumlogPiTilde = sum(logPiTilde);
  else
    d_logC0.d_alpha0 = d_logC0.d_alpha0*numgroups;
    sumlogPiTilde = cell2mat(logPiTilde);
    sumlogPiTilde = sum(sumlogPiTilde(:));    
  end
  d_Lt.d_alpha0    = d_logC0.d_alpha0 + sumlogPiTilde;
  
  %% epsilon0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  d_logC0.d_epsilon0 = K*psi(0,K*epsilon0) - K*psi(0,epsilon0);
  if ~usegroups
    sumlogATilde = sum(logATilde(:));
  else
    d_logC0.d_epsilon0 = d_logC0.d_epsilon0 * numgroups;
    sumlogATilde = cell2mat(logATilde);
    sumlogATilde = sum(sumlogATilde(:));    
  end
  d_Lt.d_epsilon0    = K*d_logC0.d_epsilon0 + sumlogATilde;
  
  %% v0 (nu0) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  d_logB0.d_v0  = 0.5*logdetW0inv - (dim/2)*log(2) - 0.5*sum(psi(0, 0.5*(v0+1-[1:dim])));
  d_Lt.d_v0     = K*d_logB0.d_v0 + 0.5*sum(logLambdaTilde);
  
  %% beta0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  d_Lt.d_beta0 = 0.5*sum(dim/beta0 - dim./beta(:) - v(:).*mWm(:));
  
  %% W0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  switch(W0mode)
    case 'iid'
      myW0inv = W0inv(1,1);
      myW0    = 1/myW0inv;
      d_logB0.d_W0  = -0.5*v0 * dim*myW0inv;
      d_trW0invW.d_W0 = zeros(1,K);
      for k=1:K
        d_trW0invW.d_W0(k) = -(myW0inv^2) * trace(W(:,:,k));
      end
      d_Lt.d_W0 = K*d_logB0.d_W0 - 0.5*sum(v(:).*d_trW0invW.d_W0(:));
    
    case 'diag'
      d_Lt.d_W0  = zeros(dim,1);
      myW0inv = diag(W0inv);
      myW0    = 1./myW0inv;
      d_logB0.d_W0  = -0.5 * v0 * myW0inv;
      d_trW0invW.d_W0 = zeros(dim,K);
      for k=1:K
        d_trW0invW.d_W0(:,k) = -(myW0inv.^2) .*diag(W(:,:,k));
      end
      d_Lt.d_W0(:) = K*d_logB0.d_W0 - 0.5*sum(repmat(v(:)',dim,1) .* d_trW0invW.d_W0, 2);
      
    otherwise
      error('bad W0mode');
  end
  
  %% m0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  tmp = zeros(dim,1);
  for k=1:K
    tmp = tmp + beta0*v(k)*W(:,:,k)*(m(:,k)-m0);
  end
  d_Lt.d_m0 = tmp;
  
  %% set derivatives to 0 when moving beyond extreme values of hyps will increase LL:
  % 1) if at max hyp, prevent increasing hyp if LL will increase
  % 2) if at min hyp, prevent decreasing hyp if LL will increase
  
  hnames = fieldnames(clipped);
  for i=1:length(hnames)
    myhname = hnames{i};
    mydname = ['d_' myhname];
    
    for j=1:length(clipped.(myhname))
      if (clipped.(myhname)(j) == +1) && (d_Lt.(mydname)(j) > 0)
        d_Lt.(mydname)(j) = 0;
      end      
      if (clipped.(myhname)(j) == -1) && (d_Lt.(mydname)(j) < 0)
        d_Lt.(mydname)(j) = 0;
      end
    end
  end
  
  %{
  if (clipped.alpha0 == +1) && (d_Lt.d_alpha0 > 0)
    d_Lt.d_alpha0 = zeros(size(d_Lt.d_alpha0));
  end  
  if (clipped.epsilon0 == +1) && (d_Lt.d_epsilon0 > 0)
    d_Lt.d_epsilon0 = zeros(size(d_Lt.d_epsilon0));
  end
  if (clipped.v0 == +1) && (d_Lt.d_v0 > 0)
    d_Lt.d_v0 = zeros(size(d_Lt.d_v0));
  end
  if (clipped.beta0 == +1) && (d_Lt.d_beta0 > 0)
    % clipping doesn't work properly w/ minimize_new?
    d_Lt.d_beta0 = zeros(size(d_Lt.d_beta0));
  end
  
  if (clipped.alpha0 == -1) && (d_Lt.d_alpha0 < 0)
    d_Lt.d_alpha0 = zeros(size(d_Lt.d_alpha0));
  end
  if (clipped.epsilon0 == -1) && (d_Lt.d_epsilon0 < 0)
    d_Lt.d_epsilon0 = zeros(size(d_Lt.d_epsilon0));
  end
  if (clipped.v0 == -1) && (d_Lt.d_v0 < 0)
    d_Lt.d_v0 = zeros(size(d_Lt.d_v0));
  end  
  if (clipped.beta0 == -1) && (d_Lt.d_beta0 < 0)
    % clipping doesn't work properly w/ minimize_new
    d_Lt.d_beta0 = zeros(size(d_Lt.d_beta0));
  end  
  
  % clip W0
  for i=1:length(d_Lt.d_W0)    
    if (clipped.W0(i) == +1) && (d_Lt.d_W0(i) > 0)
      % clipping doesn't work properly w/ minimize_new
      d_Lt.d_W0(i) = 0;
    end
    if (clipped.W0(i) == -1) && (d_Lt.d_W0(i) < 0)
      % clipping doesn't work properly w/ minimize_new
      d_Lt.d_W0(i) = 0;
    end
  end
  %}
  
  %% calculate derivatives of log (transformed) hyps
  d_LB.d_logalpha0   = d_Lt.d_alpha0   * alpha0;
  d_LB.d_logepsilon0 = d_Lt.d_epsilon0 * epsilon0;
  d_LB.d_logv0D1     = d_Lt.d_v0 * (v0-dim+1);
  d_LB.d_sqrtv0D1    = d_Lt.d_v0 * 2 * sqrt(v0-dim+1);
  d_LB.d_logbeta0    = d_Lt.d_beta0 * beta0;
  d_LB.d_sqrtbeta0   = d_Lt.d_beta0 * 2 * sqrt(beta0);
  d_LB.d_sqrtW0inv   = bsxfun(@times, d_Lt.d_W0, (myW0.^1.5) * (-2));
  d_LB.d_logW0       = bsxfun(@times, d_Lt.d_W0, myW0);
  d_LB.d_m0          = d_Lt.d_m0;
    
  % debugging
  %d_LB.LBterms = [Lt1 Lt2 Lt3 Lt4 Lt5 Lt6 Lt7 Lt8];

  
else
  d_LB = [];
  
end

  