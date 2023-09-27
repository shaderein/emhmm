function [LL_elbo,...
          sum_nu_1, ...
          update_emit_pr, ...
          update_emit_mu, ...
          update_emit_Mu,  ...
          sum_xi, ...
          sum_t_nu] = vbhem_hmm_bwd_fwd_fast(hmm_b,T,m,...
                W,v,lambda,logLambdaTilde,logATilde,logPiTilde, h3m_b_hmm_emit_cache)
 
% run the Hierarchical bwd and fwd recursions. This will compute:
% phi is the first term in the log(phi)
% theta is the E[x^{i,j}_{beta_t,rho_t-1,rho_t}] 
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2021-07-29
% Antoni B. Chan, Janet H. Hsiao, Lan Hui
% City University of Hong Kong, University of Hong Kong

% 2021-07-29: v0.80 - initial version 

g3m_b = hmm_b.emit;
M = g3m_b{1}.ncentres;

if (M>1)
  error('NOT SUPPORT for ncentres>1');
end 

% number of states 
Sb = size(hmm_b.A,1);  %row
Sr = size(m,2);
% dimension of the emission variable
dim = g3m_b{1}.nin;
   
gmmBcentres = h3m_b_hmm_emit_cache.gmmBcentres; % 1 x dim x Sb
gmmBcovars  = h3m_b_hmm_emit_cache.gmmBcovars;  % 1 x dim x Sb[diag]

%% old code 

if 0 
    tic
%get pars from hmm_b
mubl = zeros(dim,Sb);
Mubl = zeros(Sb,dim);
for i = 1:Sb
    mubl(:,i) =gmmBcentres(:,:,i)';   %column
    Mubl(i,:) = gmmBcovars(:,:,i);
end                                                    
const_denominator = (dim*log(2*pi))/2;
EElogN = zeros(Sb,Sr);


for i =1:Sb
    for j = 1:Sr       
        diff = mubl(:,i)- m(:,j);
        EElogN(i,j) = 0.5*logLambdaTilde(j) - const_denominator - 0.5*dim/lambda(j) ...
        -0.5*v(j) *( trace(W(:,:,j)*diag(Mubl(i,:))) + diff'*W(:,:,j)*diff );     
    end
end
xEElogN = EElogN;
toc
end

%% new code
    % compute the expected log-likelihood between the Gaussian components
    % i.e., E_M(r),j,rho E_M(b),i,beta E_y|M(b),I,beta [log p(y | M(r),j,rho)]
    
sum_w_pr = cell(1,Sr);
sum_w_mu = cell(1,Sr);
sum_w_Mu = cell(1,Sr);

E3logN = zeros(Sb,Sr);

for rho = 1 : Sr
    
    xfoo_sum_w_pr = zeros(Sb,M);
    xfoo_sum_w_mu = zeros(Sb,dim,M);
    switch(g3m_b{1}.covar_type)
        case 'diag'
            xfoo_sum_w_Mu = zeros(Sb,dim,M);
        case 'full'
            xfoo_sum_w_Mu = zeros(Sb,dim,dim,M);
    end
    
    
    switch(g3m_b{1}.covar_type)
        case 'diag'
            
            % sum(gmmB.covars(a,:)./gmmR.covars(b,:),2)
            ELLsall = sum(...
                bsxfun(@times, reshape(gmmBcovars,[1, dim, Sb]), ...
                reshape(repmat(diag(W(:,:,rho)),[1,Sb]), [1, dim, Sb])), 2);
            
            % + sum(((gmmB.centres(a,:) - gmmR.centres(b,:)).^2) ./ gmmR.covars(b,:),2)
            ELLsall = ELLsall + ...
                sum(bsxfun(@times, ...
                bsxfun(@minus, reshape(gmmBcentres, [ 1, dim, Sb]), ...
                reshape(m(:,rho), [ 1, dim, 1])).^2, ...
                reshape(diag(W(:,:,rho)), [ 1, dim, 1])), ...
                2);
            
            %  + dim*log(2*pi) + sum( log(gmmR.covars(b,:)),2) ...
            ELLsall =  v(rho) * ELLsall -logLambdaTilde(rho) ...
                + dim/lambda(rho) + dim*log(2*pi);
            
            % * -0.5
            ELLsall = ELLsall * -0.5;
            
        case 'full'
     
            
            
%             ELLsall = zeros(1,Sb);
%             
%             for beta = 1:Sb
%                 ELLsall_tmp = 0;
%                 ELLsall_tmp = ELLsall_tmp + trace(W(:,:,rho)*gmmBcovars(:,:,1,beta));
%                 tmp = gmmBcentres(:,:,beta) -m(:,rho)';
%                 ELLsall_tmp = ELLsall_tmp + tmp*W(:,:,rho)*tmp';
%                 ELLsall_tmp = v(rho)*ELLsall_tmp;
%                 ELLsall_tmp = ELLsall_tmp - logLambdaTilde(rho) + dim/lambda(rho)+ dim*log(2*pi);
%                 ELLsall(beta) = ELLsall_tmp;
%             end
%             ELLsall = ELLsall * -0.5;

            
            
            tmp_W = W(:,:,rho);
            trcov = reshape(sum(sum(bsxfun(@times, tmp_W, gmmBcovars), 1), 2), [1  Sb]);
            
            diff = bsxfun(@minus, gmmBcentres, m(:,rho)'); % [BK dim N]
            tmp = reshape(permute(diff, [2 1 3]), [dim 1*Sb]);
            dterm = sum(tmp .* (tmp_W * tmp), 1);
            %dterm = reshape(dterm, [1 1 1 Sb]);
            ELLsall = -0.5*( dim*log(2*pi) - logLambdaTilde(rho) + dim/lambda(rho) + ...
                v(rho)*(trcov + dterm));
                     
            
        otherwise
            error('not supported')
    end
    E3logN(:,rho) = reshape(ELLsall,[Sb,1]);
    
    
    
    xfoo_sum_w_pr = ones(Sb,1);
    
    
    xfoo_sum_w_mu = reshape(gmmBcentres,[dim,Sb])';
    
     switch(g3m_b{1}.covar_type)
       case 'diag'
         tmp3 = reshape((gmmBcentres.^2 + gmmBcovars), [1, 1, dim, Sb]);
         xfoo_sum_w_Mu = permute(reshape(sum(tmp3, 1), [1, dim, Sb]), [3 2 1]);
        
       case 'full'   
           
        cc = bsxfun(@times, reshape(gmmBcentres, [1 1 dim 1 Sb]), ...
             reshape(gmmBcentres, [1 1 1 dim Sb])); % [BK RK dim dim N]
        cc = bsxfun(@plus, cc, reshape(permute(gmmBcovars, [3 1 2 4]), [1 1 dim dim Sb]));
       
        xfoo_sum_w_Mu = reshape(sum(cc, 1), [1 dim dim Sb]);
        xfoo_sum_w_Mu = permute(xfoo_sum_w_Mu, [4 2 3 1]); 
%          
     end
%     
    sum_w_pr{rho} = xfoo_sum_w_pr;
    sum_w_mu{rho} = xfoo_sum_w_mu;
    sum_w_Mu{rho} = xfoo_sum_w_Mu; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% do the backward recursion %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
Ab = hmm_b.A;
pib= hmm_b.prior;
Ar = logATilde;
 
% allocate Theta, i.e., the VA parameters in the form of a CMC
Theta = zeros(Sr,Sr,Sb,T); % the dimensions refer to rho sigma (states of M(r)) gamma (state of M(b)) and time
 
 
LL_old  = zeros(Sb,Sr);
 
for t = T : -1 : 2
 
    LL_new = zeros(size(LL_old));
    
    if 0
    
     for rho = 1 : Sr
         
        % [N2xN] =          [1xN2]' x [1xN]      + [N,N2]'   +  [N,N2]'
        logtheta = (Ar(rho,:)' * ones(1,Sb)) + E3logN' + LL_old';  %%%phi in vhem
        % % % logtheta = (log(Ar(rho,:))' * ones(1,N)) + LLG_elbo + LL_old;
        % indexed by rho_t*beta_t
        logsumtheta = logtrick(logtheta);
 
        % normalize so that each clmn sums to 1 (may be not necessary ...) 
        theta = exp(logtheta - ones(Sr,1) * logsumtheta);  %(rho_t,beta_t)
        
        LL_new(:,rho) =  Ab * logsumtheta';
        
        % and store for later  
        Theta(rho,:,:,t) = theta;  %%\rho means \rho_{t-1}
     end
     
     xTheta = Theta;
     xLL_new = LL_new;
   end
    
      %% 2017-08-07 - ABC - fast code
  % [N2xN]xN2 =          [1xN2]' x [1xN]      + [N,N2]'   +  [N,N2]'
  %    logtheta = (Ar(rho,:)' * ones(1,Sb)) + EElogN' + LL_old'; 
    logtheta_all = bsxfun(@plus, ...
        reshape(Ar', [Sr 1 Sr]), ...
        reshape(E3logN',[Sr Sb]) + reshape(LL_old',[Sr Sb]));
 
  % [1xSb]xSr
  %logsumtheta = logtrick(logtheta);
    logsumtheta_all = logtrick2(logtheta_all, 1);
    
  % [NxSr]         =  [NxN] * [1xN]'
  %  LL_new(:,rho) =  Ab * logsumtheta';
    LL_new = reshape(sum(bsxfun(@times, reshape(Ab,[Sb Sb 1]), logsumtheta_all), 2), [Sb Sr]);
  
  % normalize so that each clmn sums to 1 (may be not necessary ...)
  % [N2xN]xN2 = [N2xN]   -       [N2x1]*[1xN]
  % theta = exp(logtheta - ones(Sr,1) * logsumtheta);
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
 
logtheta = (logPiTilde * ones(1,Sb)) + E3logN' + LL_old';
 
logsumtheta = logtrick(logtheta);
 
LL_elbo =  pib' * logsumtheta';
 
theta = exp(logtheta - ones(Sr,1) * logsumtheta); 
 
Theta_1 = theta;

%%%%%%%%%%%%%%  add new %%%%
% Theta_all.Theta_1 = Theta_1;
% Theta_all.Theta = Theta;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% do the forward recursion  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% rather then saving all intermediate values, just cache the cumulative
% values that are needed for the updates (it saves a lot of memory)
 
% nu is N2 by N, first dimension indexed by sigma (M(r)), second by gamma (M(b))
 
% initialize
nu = (ones(Sr,1) * pib') .* Theta_1; % N2 by N (indexed by sigma and gamma)
%xnu=nu;
sum_nu_1 = sum(nu,2)';    %hat v_1
sum_t_nu = nu;
% CACHE: sum_t sum_gamma xi(rho,sigma,gamma,t)
%xsum_t_sum_g_xi = zeros(Sr,Sr); % N2 by N2 (indexed by rho and sigma) hat_xi
sum_t_sum_g_xi = zeros(Sr,Sr); 
xsum_t_sum_g_xi = zeros(Sr,Sr);
%%try 3.12
%nu_theta_1 = sum(sum(nu.*log(Theta_1)));  %try
%sum_t_nu_theta = nu_theta_1;  %try
     
% CACHE: sum_t sum_gamma xi(rho,sigma,gamma,t)
 
for t = 2 : T
    % compute the inner part of the update of xi (does not depend on sigma)
    foo = nu * Ab; % indexed by rho gamma
 
    if 0
    for sigma = 1 : Sr
 
        % ABC: bug fix when another dim is 1
        temp_T = reshape(Theta(:,sigma,:,t), [size(Theta,1), size(Theta,3)]);
        xi_foo = foo .*temp_T; % (indexed by rho gamma);
        
        %%try
        %xi_foo_log_phi = xi_foo.*log(reshape(Theta(:,sigma,:,t), [size(Theta,1), size(Theta,3)])+1e-50);
        %sum_t_nu_theta =sum_t_nu_theta  + sum(sum(xi_foo_log_phi));
        
        % CACHE:
        xsum_t_sum_g_xi(:,sigma) = xsum_t_sum_g_xi(:,sigma) + sum(xi_foo,2);       
        % new nu
        % nu(sigma,:) = ones(1,N) * squeeze(xi(:,sigma,:,t));
        xnu(sigma,:) = ones(1,Sr) * xi_foo;      
    end
    end
 
 
 
 if 1
    %% 2018-08-07 - ABC - faster code w/o loop
    % xi_foo = foo .* reshape(Theta(:,sigma,:,t), [size(Theta,1), size(Theta,3)]); % (indexed by rho gamma);
    xi_foo_all = bsxfun(@times, reshape(foo,[Sr 1 Sb]), Theta(:,:,:,t));
    
    % sum_t_sum_g_xi(:,sigma) = sum_t_sum_g_xi(:,sigma) + sum(xi_foo,2);
    sum_t_sum_g_xi = sum_t_sum_g_xi + sum(xi_foo_all, 3);
    
    % nu(sigma,:) = ones(1,N2) * xi_foo;
    nu = reshape(sum(xi_foo_all,1), [Sr Sb]);
 end
    
    %% consistency check
    if 0
      err = [];
      err(end+1) = sum(sum(abs(xnu-nu)));
      err(end+1) = sum(sum(abs(xsum_t_sum_g_xi - sum_t_sum_g_xi)));
      if any(err>1e-10)
        warning('mismatch error');
        keyboard
      end
    end 
                           
    % CACHE: in the sum_t nu_t(sigma, gamma)
    sum_t_nu = sum_t_nu + nu;  
     
end
 
%LL_elbo = LL_elbo/T;
sum_xi = sum_t_sum_g_xi;
%sum_theta =sum_t_nu_theta;

 

%%
%%%% now prepare the cumulative sufficient statistics for the reestimation
%%%% of the emission distributions

update_emit_pr = zeros(Sr,M);
update_emit_mu = zeros(Sr,dim,M);
switch hmm_b.emit{1}.covar_type
  case 'diag'
    update_emit_Mu = zeros(Sr,dim,M);
  case 'full'
    update_emit_Mu = zeros(Sr,dim,dim,M);
end


% loop all the emission GMM of each state
for sigma = 1 : Sr
  
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
          reshape(sum_t_nu(sigma,:), [Sb 1 1]), foo_sum_w_Mu(:,:,:,l)), 1);
    end
  end
  
  
end
% %% symmary the output
%  
%  output.sum_nu_1 = sum_nu_1;
%  output.sum_xi = sum_xi;
%  %output.sum_theta = sum_theta;
%  output.sum_t_nu = sum_t_nu;
%  output.Theta = Theta;
%  output.update_emit_pr = update_emit_pr;
%  output.update_emit_mu = update_emit_mu;
%  output.update_emit_Mu = update_emit_Mu;
 
 %output.Theta_all = Theta_all;
 
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


 

