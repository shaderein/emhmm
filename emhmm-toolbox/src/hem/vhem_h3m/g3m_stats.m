function [LLG_elbo, sum_w_pr, sum_w_mu, sum_w_Mu, sum_w_mu2] = g3m_stats(g3m_b,g3m_r, h3m_b_hmm_emit_cache)
% H3M Toolbox 
% Copyright 2012, Emanuele Coviello, CALab, UCSD

% 2017-08-06 - ABC - MATLAB speed optimizations
% 2017-08-17 - ABC - added full covariance matrices
% 2018-05-25: v0.73 - added "sum_w_mu2" output for MANOVA test.

% first, get the elbo of the E_gauss(b) [log p (y | gauss (r))]  (different b on different rows)
% and the sufficient statistics for the later updates
% sum_w_pr is a N by 1 cell of N by M
% sum_w_mu is a N by 1 cell of N by dim by M
% sum_w_Mu is a N by 1 cell of N by dim by M (for the diagonal case)
% sum_w_Mu is a N by 1 cell of N by dim by dim by M (for the full case)
%
% if necessary:
% sum_w_mu2 is a N by 1 cell of N by dim by dim by M

% v0.73 - compute mu2?
if nargout<5
  COMPUTE_MU2 = 0;
else
  COMPUTE_MU2 = 1;
end

% number of states, e.g., how many GMM emission pdf
N = length(g3m_b);
N2 = length(g3m_r);
M = g3m_b{1}.ncentres;
dim = g3m_b{1}.nin;

if (M>1)
  error('optimized g3m_stats has not been tested for ncentres>1');
end
  
%% old code %%
if 0
  
  LLG_elbo = zeros (N,N2);
  
  % compute variational elbo to E_M(b),beta [log p(y | M(r),rho)]
  % and variational posteriors
  sum_w_pr = cell(1,N2);
  sum_w_mu = cell(1,N2);
  sum_w_Mu = cell(1,N2);
  
  for rho = 1 : N2
    
    
    foo_sum_w_pr = zeros(N,M);
    foo_sum_w_mu = zeros(N,dim,M);
    switch(g3m_b{1}.covar_type)
      case 'diag'    
        foo_sum_w_Mu = zeros(N,dim,M);
      case 'full'
        foo_sum_w_Mu = zeros(N,dim,dim,M);
    end
    
    gmmR = g3m_r{rho};
    
    for beta = 1 : N
      
      gmmB = g3m_b{beta};
      
      % compute the expected log-likelihood between the Gaussian components
      % i.e., E_M(b),beta,m [log p(y | M(r),rho,l)], for m and l 1 ...M
      [ELLs] =  compute_exp_lls(gmmB,gmmR);
      
      % compute log(omega_r) + E_M(b),beta,m [log p(y | M(r),rho,l)]
      log_theta = ELLs + ones(M,1) * log(gmmR.priors);
      
      % compute log Sum_b omega_b exp(-D(fa,gb))
      log_sum_theta = logtrick(log_theta')';
      
      % compute L_variational(M(b)_i,M(r)_j) = Sum_a pi_a [  log (Sum_b omega_b exp(-D(fa,gb)))]
      
      LLG_elbo(beta,rho) = gmmB.priors * log_sum_theta;
      
      % cache theta
      theta = exp(log_theta -  log_sum_theta * ones(1,M));
      %          Theta{beta,rho} = theta; % indexed by m and l, [M by M]
      
      % aggregate in the output ss
      
      foo_sum_w_pr(beta,:) = gmmB.priors * theta;
      
      %          for l = 1 :  M
      %
      %              foo_sum_w_mu(beta,:,l) = foo_sum_w_mu(beta,:,l) + (gmmB.priors .* theta(:,l)') * gmmB.centres;
      %
      %          end
      
      foo_sum_w_mu(beta,:,:) =  ( ((ones(M,1) * gmmB.priors) .* theta') * gmmB.centres )';
      
      switch(g3m_b{1}.covar_type)
        case 'diag'
          foo_sum_w_Mu(beta,:,:) =  ( ((ones(M,1) * gmmB.priors) .* theta') * (gmmB.centres.^2 + gmmB.covars) )';
        
        case 'full'
          tmp = gmmB.covars;
          for m=1:M
            tmp(:,:,m) = tmp(:,:,m) + gmmB.centres(m,:)'*gmmB.centres(m,:);
          end
          
          ww = ((ones(M,1) * gmmB.priors) .* theta'); %[M x M]
          
          for m=1:M
            foo_sum_w_Mu(beta,:,:,m) = sum(bsxfun(@times, tmp, reshape(ww(m,:), [1 1 M])), 3);
          end
      end
      
    end
    
    sum_w_pr{rho} = foo_sum_w_pr;
    sum_w_mu{rho} = foo_sum_w_mu;
    sum_w_Mu{rho} = foo_sum_w_Mu;
  end
  unopt_time = toc;
  
  % save for comparison
  x_sum_w_pr = sum_w_pr;
  x_sum_w_mu = sum_w_mu;
  x_sum_w_Mu = sum_w_Mu;
  x_LLG_elbo = LLG_elbo;
  
end
  


if 1
  
  %% ABC - 2017-08-06 - faster code without inner loop %%
  LLG_elbo = zeros (N,N2);
  
  % compute variational elbo to E_M(b),beta [log p(y | M(r),rho)]
  % and variational posteriors
  sum_w_pr = cell(1,N2);
  sum_w_mu = cell(1,N2);
  sum_w_Mu = cell(1,N2);  
  if COMPUTE_MU2
    sum_w_mu2 = cell(1,N2); % v0.73
  end
  
  dim = g3m_b{1}.nin;
  BKmax = g3m_b{1}.ncentres;
  
  gmmBcentres = h3m_b_hmm_emit_cache.gmmBcentres; % BKmax x dim x N
  gmmBcovars  = h3m_b_hmm_emit_cache.gmmBcovars;  % BKmax x dim x N [diag]
                                                  % dim x dim x BKmax x N [full]
  gmmBpriors  = h3m_b_hmm_emit_cache.gmmBpriors;  % BKmax x N
  
  
  
  for rho = 1 : N2
    
    xfoo_sum_w_pr = zeros(N,M);
    xfoo_sum_w_mu = zeros(N,dim,M);
    switch(g3m_b{1}.covar_type)
      case 'diag'
        xfoo_sum_w_Mu = zeros(N,dim,M);
      case 'full'
        xfoo_sum_w_Mu = zeros(N,dim,dim,M);
    end    
    if COMPUTE_MU2
      xfoo_sum_w_mu2 = zeros(N,dim,dim,M); % v0.73
    end
    
    gmmR = g3m_r{rho};
    
    
    RK = gmmR.ncentres;
    
    % compute the expected log-likelihood between the Gaussian components
    % i.e., E_M(b),beta,m [log p(y | M(r),rho,l)], for m and l 1 ...M
    %[ELLs] =  compute_exp_lls(gmmB,gmmR);
    
    
    % ELLsall = [BKmax x RK x 1 x N]
    
    switch(g3m_b{1}.covar_type)
      case 'diag'
        %
        %    ELLs = zeros(A,B);
        %    for a=1:A
        %      for b=1:B
        %        ELLs(a,b) = -.5 * ( ...
        %                     dim*log(2*pi) + sum( log(gmmR.covars(b,:)),2) ...
        %                     +  sum(gmmB.covars(a,:)./gmmR.covars(b,:),2)  ...
        %                     + sum(((gmmB.centres(a,:) - gmmR.centres(b,:)).^2) ./ gmmR.covars(b,:),2) ...
        %                    );
        %      end
        %    end
        
        sumlogRcovars = sum( log(gmmR.covars),2);
        
        % sum(gmmB.covars(a,:)./gmmR.covars(b,:),2)
        ELLsall = sum(...
          bsxfun(@rdivide, reshape(gmmBcovars,[BKmax, 1, dim, N]), ...
          reshape(gmmR.covars, [1, RK, dim, 1])), 3);
        
        % + sum(((gmmB.centres(a,:) - gmmR.centres(b,:)).^2) ./ gmmR.covars(b,:),2)
        ELLsall = ELLsall + ...
          sum(bsxfun(@rdivide, ...
          bsxfun(@minus, reshape(gmmBcentres, [BKmax, 1, dim, N]), ...
          reshape(gmmR.centres, [1, RK, dim, 1])).^2, ...
          reshape(gmmR.covars, [1, RK, dim, 1])), ...
          3);
        
        %  + dim*log(2*pi) + sum( log(gmmR.covars(b,:)),2) ...
        ELLsall = bsxfun(@plus, ELLsall, reshape(sumlogRcovars, [1 RK 1 1])) ...
          + dim*log(2*pi);
        
        % * -0.5
        ELLsall = ELLsall * -0.5;
        
      case 'full'        
        % ELLs = zeros(A,B);
        % for a=1:A
        %   for b=1:B
        %     inv_covR = inv(gmmR.covars(:,:,b));
        %    
        %     ELLs(a,b) = -.5 * ( ...
        %       dim*log(2*pi) + log(det(gmmR.covars(:,:,b))) ...
        %        + trace( inv_covR * gmmB.covars(:,:,a) ) +  ...
        %        + (gmmB.centres(a,:) - gmmR.centres(b,:)) * inv_covR * (gmmB.centres(a,:) - gmmR.centres(b,:))' ...
        %      );
        %   end
        % end
        
        ELLsall = zeros(BKmax, RK, 1, N);
        for j=1:RK
          inv_covR = inv(gmmR.covars(:,:,j));
          trcov = reshape(sum(sum(bsxfun(@times, inv_covR, gmmBcovars), 1), 2), [BKmax 1 1 N]);
                    
          diff = bsxfun(@minus, gmmBcentres, gmmR.centres(j,:)); % [BK dim N]
          tmp = reshape(permute(diff, [2 1 3]), [dim BKmax*N]);
          dterm = sum(tmp .* (inv_covR * tmp), 1);
          dterm = reshape(dterm, [BKmax 1 1 N]);
          ELLsall(:,j,1,:) = -0.5*( dim*log(2*pi) + log(det(gmmR.covars(:,:,j))) + ...
            trcov + dterm);
        end        
        
        if 0
          % consistency check
          ELLsx = zeros(BKmax, RK, 1, N);
          for beta=1:N
            ELLsx(:,:,1,beta) =  compute_exp_lls(g3m_b{beta},gmmR);
          end
          err = sum(abs(ELLsall(:) - ELLsx(:)));
          if (err>1e-6)
            keyboard
          end
        end
    end
    
    
    % compute log(omega_r) + E_M(b),beta,m [log p(y | M(r),rho,l)]
    % log_theta = ELLs + ones(M,1) * log(gmmR.priors);
    %  [BKmax x RK x 1 x N]
    log_theta_all = bsxfun(@plus, ELLsall, reshape(log(gmmR.priors), [1,RK,1,1]));
    
    % compute log Sum_b omega_b exp(-D(fa,gb))
    % log_sum_theta = logtrick(log_theta')';
    %  [BKmax x 1 x 1 x N]
    log_sum_theta_all = logtrick2(log_theta_all, 2);
    
    % compute L_variational(M(b)_i,M(r)_j) = Sum_a pi_a [  log (Sum_b omega_b exp(-D(fa,gb)))]
    % LLG_elbo(beta,rho) = gmmB.priors * log_sum_theta;
    LLG_elbo(:,rho) = sum(bsxfun(@times, reshape(gmmBpriors,[BKmax, N]), ...
      reshape(log_sum_theta_all, [BKmax, N])), 1);
    
    % cache theta
    % theta = exp(log_theta -  log_sum_theta * ones(1,M));
    % [BKmax x RK x 1 x N]
    theta_all = exp(bsxfun(@minus, log_theta_all, log_sum_theta_all));
    
    % aggregate in the output ss
    % foo_sum_w_pr(beta,:) = gmmB.priors * theta;
    % [N x RK]
    xfoo_sum_w_pr = reshape(sum(...
      bsxfun(@times, reshape(gmmBpriors, [BKmax, 1, 1, N]), ...
      theta_all), 1), [RK N])';
    
    %xfoo_sum_w_mu(beta,:,:) =  ( ((ones(M,1) * gmmB.priors) .* theta') * gmmB.centres )';
    tmp  = bsxfun(@times, reshape(gmmBpriors, [BKmax, 1, 1, N]), theta_all);
    tmp2 = bsxfun(@times, tmp, reshape(gmmBcentres, [BKmax, 1, dim, N]));
    xfoo_sum_w_mu = permute(reshape(sum(tmp2, 1), [RK, dim, N]), [3 2 1]);
    
    %xfoo_sum_w_Mu(beta,:,:) =  ( ((ones(M,1) * gmmB.priors) .* theta') * (gmmB.centres.^2 + gmmB.covars) )';
    switch(g3m_b{1}.covar_type)
      case 'diag'
        tmp3 = bsxfun(@times, tmp, reshape((gmmBcentres.^2 + gmmBcovars), [BKmax, 1, dim, N]));
        xfoo_sum_w_Mu = permute(reshape(sum(tmp3, 1), [RK, dim, N]), [3 2 1]);
        
        
      case 'full'                
        cc = bsxfun(@times, reshape(gmmBcentres, [BKmax 1 dim 1 N]), ...
             reshape(gmmBcentres, [BKmax 1 1 dim N])); % [BK RK dim dim N]
        cc = bsxfun(@plus, cc, reshape(permute(gmmBcovars, [3 1 2 4]), [BKmax 1 dim dim N]));
        tmp3 = bsxfun(@times, reshape(tmp, [BKmax RK 1 1 N]), cc);
        xfoo_sum_w_Mu = reshape(sum(tmp3, 1), [RK dim dim N]);
        xfoo_sum_w_Mu = permute(xfoo_sum_w_Mu, [4 2 3 1]); 
        % [zeros(N,dim,dim,M)]
         
    end
    
    if COMPUTE_MU2
      % v0.73: xfoo_sum_w_mu2(beta,:,:) =  ( ((ones(M,1) * gmmB.priors) .* theta') * (gmmB.centres.^2)
      cc = bsxfun(@times, reshape(gmmBcentres, [BKmax 1 dim 1 N]), ...
        reshape(gmmBcentres, [BKmax 1 1 dim N])); % [BK RK dim dim N]
      tmp3 = bsxfun(@times, reshape(tmp, [BKmax RK 1 1 N]), cc);
      xfoo_sum_w_mu2 = reshape(sum(tmp3, 1), [RK dim dim N]);
      xfoo_sum_w_mu2 = permute(xfoo_sum_w_mu2, [4 2 3 1]);
      sum_w_mu2{rho} = xfoo_sum_w_mu2; % v0.73
    end
    
    sum_w_pr{rho} = xfoo_sum_w_pr;
    sum_w_mu{rho} = xfoo_sum_w_mu;
    sum_w_Mu{rho} = xfoo_sum_w_Mu;          
  end
end

if 0
  % consistency check
  err(1) = sum(abs(LLG_elbo(:) - x_LLG_elbo(:)));
  err(2) = sum(sum(sum(sum(abs(cat(1,sum_w_pr{:}) - cat(1,x_sum_w_pr{:}))))));
  err(3) = sum(sum(sum(sum(abs(cat(1,sum_w_mu{:}) - cat(1,x_sum_w_mu{:}))))));
  err(4) = sum(sum(sum(sum(abs(cat(1,sum_w_Mu{:}) - cat(1,x_sum_w_Mu{:}))))));
  
  if any(err>1e-6)
    keyboard
  end
end


end

function [ELLs] = compute_exp_lls(gmmA,gmmB)

dim = gmmA.nin;
A = gmmA.ncentres;
B = gmmB.ncentres;

ELLs = zeros(A,B);

for a = 1 : A
  for b = 1 : B
    
    switch(gmmA.covar_type)
      %                 case 'spherical'
      %                     KLs(a,b) = .5 * ( ...
      %                          dim*log( gmmB.covars(1,b) / gmmA.covars(1,a) )...
      %                          +  dim* gmmA.covars(1,a)/gmmB.covars(1,b) - dim ...
      %                          + sum(((gmmA.centres(a,:) - gmmB.centres(b,:)).^2) ./ gmmB.covars(1,b)) ...
      %                         );
      
      case 'diag'
        %                     KLs(a,b) = .5 * ( ...
        %                          sum( log(gmmB.covars(b,:)),2) - sum(log(gmmA.covars(a,:)),2 )...
        %                          +  sum(gmmA.covars(a,:)./gmmB.covars(b,:),2) - dim ...
        %                          + sum(((gmmA.centres(a,:) - gmmB.centres(b,:)).^2) ./ gmmB.covars(b,:),2) ...
        %                         );
        ELLs(a,b) = -.5 * ( ...
          dim*log(2*pi) + sum( log(gmmB.covars(b,:)),2) ...
          +  sum(gmmA.covars(a,:)./gmmB.covars(b,:),2)  ...
          + sum(((gmmA.centres(a,:) - gmmB.centres(b,:)).^2) ./ gmmB.covars(b,:),2) ...
          );
        
      case 'full'
        inv_covB = inv(gmmB.covars(:,:,b));
        
        %KLs(a,b) = .5 * ( ...
        %  log(det(gmmB.covars(:,:,b))) - log(det(gmmA.covars(:,:,a)))...
        %  + trace( inv_covB * gmmA.covars(:,:,a) ) - dim ...
        %  + (gmmA.centres(a,:) - gmmB.centres(b,:)) * inv_covB * (gmmA.centres(a,:) - gmmB.centres(b,:))' ...
        %  );
        ELLs(a,b) = -.5 * ( ...
          dim*log(2*pi) + log(det(gmmB.covars(:,:,b))) ...
          + trace( inv_covB * gmmA.covars(:,:,a) ) +  ...
          + (gmmA.centres(a,:) - gmmB.centres(b,:)) * inv_covB * (gmmA.centres(a,:) - gmmB.centres(b,:))' ...
          );
        
      otherwise
        error('Covarance type not supported')
    end
    
  end
  
end

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