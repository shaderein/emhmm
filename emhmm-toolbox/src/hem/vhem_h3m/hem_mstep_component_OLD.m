function hmm = hem_mstep_component(hmm, estats, mopt)
% hem_mstep_component - compute M-step for one component
%
%  hmm = hem_mstep_component(hmm, estats, mopt)
% 
%  INPUTS:
%  hmm - template for hmm
%
%  E-step statistics:
%      estats.Zomega  = Z(:,j) .* h3m_b.omega(:);
%      estats.nu_1    = nu_1(:,j);           (cell)
%      estats.sum_xi  = sum_xi(:,j);         (cell)
%      estats.emit_pr = update_emit_pr(:,j); (cell) 
%      estats.emit_mu = update_emit_mu(:,j); (cell)
%      estats.emit_M  = update_emit_M(:,j);  (cell)
%
%      estats.emit_mu2 = update_emit_mu2(:,j); (cell) [optional]
%
%  constants:
%      estats.M   = M;
%      estats.N2  = size(h3m_r.hmm{j}.A,1);
%      estats.dim = dim;
%      estats.Kb  = Kb;
%      estats.reg_cov = reg_cov;
%
%  mopt - options
%
% OUTPUT
%  hmm - the updated HMM component
%        (stats will be computed if estats.emit_mu2 is present)
%
% ---
% H3M Toolbox 
% Copyright 2012, Emanuele Coviello, CALab, UCSD
% Antoni Chan 2017/06/06

% 2017-08-17: ABC: added full covariance matrices
% 2018-05-25: v0.73 - output mu2, stats structure
% 2018-06-08: v0.73 - fix bug with A when only one state
% 2018-06-13: v0.73 - fix bug with A when tau=1

% v0.73 - check to compute mu2
if isfield(estats, 'emit_mu2')
  COMPUTE_MU2 = 1;
else
  COMPUTE_MU2 = 0;
end

M   = estats.M;
N2  = estats.N2;
dim = estats.dim;
Kb  = estats.Kb;

% initialize matrices
new_prior = zeros(N2,1);
new_A     = zeros(N2,N2);
        
new_Gweight = cell(1,N2);
new_Gmu     = cell(1,N2);
new_GMu     = cell(1,N2);
        
for n = 1 : N2
  new_Gweight{n}  = zeros(1,M);
  new_Gmu{n}      = zeros(M,dim);
            
  switch mopt.emit.covar_type
    case 'diag'
      new_GMu{n}       = zeros(M,dim);
    case 'full'
      new_GMu{n}       = zeros(dim,dim,M);
  end    
end

% v0.73 - initialize mu2
if COMPUTE_MU2
  new_Gmu2 = cell(1,N2);
  for n=1:N2
    new_Gmu2{n} = zeros(dim,dim,M);
  end
end
  

% pull out commonly used terms
Zj = estats.Zomega;

% loop all the components of the base mixture
for i = 1 : Kb                       
  if Zj(i) > 0
            
    nu     = estats.nu_1{i};           % this is a 1 by N vector
    xi     = estats.sum_xi{i};         % this is a N by N matrix (from - to)
    up_pr  = estats.emit_pr{i}; % this is a N by M matrix
    up_mu  = estats.emit_mu{i}; % this is a N by dim by M matrix
    up_Mu  = estats.emit_M{i};  % this is a N by dim by M matrix [diagonal covariance]
                                % this is a N by dim by dim by M matrix [full covariance]
    
    new_prior = new_prior + Zj(i) * nu';
    new_A     = new_A     + Zj(i) * xi;

    for n = 1 : N2
      new_Gweight{n}  = new_Gweight{n} + Zj(i) * up_pr(n,:);
      % % %                 new_Gmu{n}      = new_Gmu{n}     + Z(i,j) * squeeze(up_mu(n,:,:))';
      new_Gmu{n}      = new_Gmu{n}     + Zj(i) * reshape(up_mu(n,:,:),dim,[])';

      switch mopt.emit.covar_type
        case 'diag'
          % % %                     new_GMu{n}  = new_GMu{n}     + Z(i,j) * squeeze(up_Mu(n,:,:))';
          new_GMu{n}  = new_GMu{n}     + Zj(i) * reshape(up_Mu(n,:,:),dim,[])';
          
        case 'full'
          new_GMu{n}  = new_GMu{n} + Zj(i) * reshape(up_Mu(n,:,:,:), [dim dim M]);
      end
      
      % v0.73 - accumulate mu2
      if COMPUTE_MU2
        new_Gmu2{n} = new_Gmu2{n} + Zj(i) * reshape(estats.emit_mu2{i}(n,:,:,:), [dim dim M]);
      end
      
    end
  end            
end
        
% BUG FIX: v0.73: if only 1 state, then xi=0. So force A to be a small number.
if (N2==1)
  new_A = 1e-12;
end

% BUG FIX: v0.73: if tau=1, then xi=0, so force A to be identity matrix (no transitions)
if (mopt.tau == 1)
  new_A = 1e-12*eye(size(new_A));
end

% normalize things, i.e., divide by the denominator
hmm.prior = new_prior / sum(new_prior);
hmm.A     = new_A    ./ repmat(sum(new_A,2),1,N2);
        
% ABC 2016-12-09 - save the counts in each emission
hmm.stats.emit_vcounts = sum(new_A,1) + new_prior';

% normalize the emission distrbutions
for n = 1 : N2
            
  %normalize the mean
  hmm.emit{n}.centres = new_Gmu{n} ./ (new_Gweight{n}' * ones(1,dim));
            
  %normalize the covariance
  switch mopt.emit.covar_type
    case 'diag'
      Sigma  = new_GMu{n}  - 2* (new_Gmu{n} .* hmm.emit{n}.centres) ...
        + (hmm.emit{n}.centres.^2) .* (new_Gweight{n}' * ones(1,dim));
      hmm.emit{n}.covars = Sigma ./ (new_Gweight{n}' * ones(1,dim));
      hmm.emit{n}.covars = hmm.emit{n}.covars + estats.reg_cov;
      
    case 'full'
      for m=1:M
        tmp  = new_Gmu{n}(m,:)'*hmm.emit{n}.centres(m,:);
        tmp2 = new_Gweight{n}(m) * (hmm.emit{n}.centres(m,:)'*hmm.emit{n}.centres(m,:));
        Sigma  = new_GMu{n}(:,:,m)  - tmp - tmp' + tmp2;
        hmm.emit{n}.covars(:,:,m) = Sigma / new_Gweight{n}(m);
        hmm.emit{n}.covars(:,:,m) = hmm.emit{n}.covars(:,:,m) + estats.reg_cov*eye(dim);
        
        % test diagional full matrix
        %hmm.emit{n}.covars(:,:,m) = diag(diag(hmm.emit{n}.covars(:,:,m)));
        
      end
  end
         
  % v0.73 - normalize mu2, collect stats
  if COMPUTE_MU2
    for m=1:M
      hmm.stats.emit{n}.mu2(:,:,m) =  new_Gmu2{n}(:,:,m) / new_Gweight{n}(m);
    end    
    hmm.stats.emit{n}.mu = hmm.emit{n}.centres;
  end
    
  
  % normalize the mixture weight of the gaussian
  if (M == 1)
    hmm.emit{n}.priors = 1;
  else
    hmm.emit{n}.priors = new_Gweight{n} ./ sum(new_Gweight{n});
  end
            
  %%%%%%% if some of the emission prob is zero, replace it ...  
  %ind_zero = find(hmm.emit{n}.counts_emit == 0);
  %for i_z = ind_zero
  %  %fprintf('!!! modifying gmm emission: one component is zero \n')
  %  [foo highest] = max(hmm.emit{n}.priors);
  %  hmm.emit{n}.priors([i_z highest]) = hmm.emit{n}.priors(highest)/2;
  %  % renormalize for safety
  %  hmm.emit{n}.priors = hmm.emit{n}.priors / sum(hmm.emit{n}.priors);
  %     
  %  hmm.emit{n}.centres(i_z,:) = hmm.emit{n}.centres(highest,:);
  %  hmm.emit{n}.covars(i_z,:)  = hmm.emit{n}.covars(highest,:);
  %  
  %  % perturb only the centres
  %  hmm.emit{n}.centres(i_z,:) =  hmm.emit{n}.centres(i_z,:) + (0.01 * rand(size(hmm.emit{n}.centres(i_z,:))))...
  %          .*hmm.emit{n}.centres(i_z,:);
  %end             
             
end

