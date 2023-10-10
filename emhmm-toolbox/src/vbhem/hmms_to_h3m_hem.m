function [h3m] = hmms_to_h3m_hem(hmms, covmode,use_post)
% hmms_to_h3m_hem - convert a list of HMMs into H3M toolbox format
%
%     h3m = hmms_to_h3m(hmms, covmode)
%
% INPUT:  hmms    = a cell array of hmms
%         covmode = 'diag' - use diagonal covariance
%                 = 'full' - use full covariance
% OUTPUT: h3m  = HMM mixture format for H3M toolbox
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2021-07-29
% Antoni B. Chan, Janet H. Hsiao, Lan Hui
% City University of Hong Kong, University of Hong Kong

% 2021-07-29: v0.80 - initial version


if nargin<3
    use_post = 0;
end

% get number of HMMs
h3m.K = length(hmms);

% 2018-11-23: v0.74 - do the normalization later;
h3m.omega = ones(1,h3m.K);

% input dimension
% 2018-11-23: v0.74 - handle empty hmm
for j=1:h3m.K
  if ~isempty(hmms{j})
    nin = length(hmms{j}.pdf{1}.mean);
    break;
  end
end


% convert each HMM to the H3M toolbox format
for j=1:h3m.K  
  clear tempHmm

  if ~isempty(hmms{j})    
    % state prior, transition matrix
    S = length(hmms{j}.prior);
       
    if use_post      
        alpha = hmms{j}.varpar.alpha;   
        epsilon = hmms{j}.varpar.epsilon;%+0.001;
        
        psiEpsilonHat = zeros(1,S);
        logATilde     = zeros(S,S);
        for k = 1:S
        %Bishop (10.66)
            psiEpsilonHat(k) = psi(0,sum(epsilon(k,:)));
            logATilde(k,:) = psi(0,epsilon(k,:)) - psiEpsilonHat(k);  
        end
        
        psiAlphaHat = psi(0,sum(alpha));
        logPiTilde = psi(0,alpha) - psiAlphaHat;
        
        tempHmm.prior = exp(logPiTilde);
        tempHmm.A = exp(logATilde);
%         tempHmm.prior = exp(logPiTilde)./sum(exp(logPiTilde));
%         tempHmm.A = exp(logATilde)./sum(exp(logATilde),2);
    else
        
        tempHmm.prior = hmms{j}.prior;
        tempHmm.A     = hmms{j}.trans;
         
    end
    
    % for each emission density, convert to H3M format
    for i = 1:S
      tempHmm.emit{i}.type     = 'gmm';  % a GMM with one component
      tempHmm.emit{i}.nin      = nin;
      tempHmm.emit{i}.ncentres = 1;
      tempHmm.emit{i}.priors   = 1;
      tempHmm.emit{i}.centres  = hmms{j}.pdf{i}.mean;
      
      if use_post
          
          switch(covmode)  
              
            case 'diag'
              tempHmm.emit{i}.covar_type = 'diag';
              tilde_beta = (hmms{j}.varpar.beta(i)+1)/hmms{j}.varpar.beta(i);
              tempHmm.emit{i}.covars = tilde_beta*diag(hmms{j}.pdf{i}.cov)';
              %tempHmm.emit{i}.covars(1,1) = hmms{j}.pdf{i}.cov(1,1);
              %tempHmm.emit{i}.covars(1,2) = hmms{j}.pdf{i}.cov(2,2);
            case 'full'
              tempHmm.emit{i}.covar_type = 'full';
              tilde_beta = (hmms{j}.varpar.beta(i)+1)/hmms{j}.varpar.beta(i);
              tempHmm.emit{i}.covars     = tilde_beta.*hmms{j}.pdf{i}.cov;
            otherwise
              error('unknown covmode')
          end
          
      else          
           
          switch(covmode)
            case 'diag'
              tempHmm.emit{i}.covar_type = 'diag';
              tempHmm.emit{i}.covars = diag(hmms{j}.pdf{i}.cov)';
              %tempHmm.emit{i}.covars(1,1) = hmms{j}.pdf{i}.cov(1,1);
              %tempHmm.emit{i}.covars(1,2) = hmms{j}.pdf{i}.cov(2,2);
            case 'full'
              tempHmm.emit{i}.covar_type = 'full';
              tempHmm.emit{i}.covars     = hmms{j}.pdf{i}.cov;

              %tempHmm.emit{i}.covars     = diag(diag(hmms{j}.pdf{i}.cov)); % test diagonal full matrix
            otherwise
              error('unknown covmode')
          end      
      end     
    end
    
  else
    % 2018-11-23: v0.74 - make a dummy HMM w/ one state, and make the weight 0
    tempHmm.prior = 1;
    tempHmm.A     = 1;    
    tempHmm.emit{1}.type     = 'gmm';  % a GMM with one component
    tempHmm.emit{1}.nin      = nin;
    tempHmm.emit{1}.ncentres = 1;
    tempHmm.emit{1}.priors   = 1;
    tempHmm.emit{1}.centres  = zeros(1,nin);    
    switch(covmode)
      case 'diag'
        tempHmm.emit{1}.covar_type = 'diag';
        tempHmm.emit{1}.covars = ones(1,nin);
        case 'full'
        tempHmm.emit{1}.covar_type = 'full';
        tempHmm.emit{1}.covars     = eye(nin);
      otherwise
        error('unknown covmode')
    end
    h3m.omega(j) = 0;
  end
    
  % save to new structure
  h3m.hmm{j} = tempHmm;
end

% 2018-11-23: v0.74 - now normalize the weights
h3m.omega = h3m.omega / sum(h3m.omega);


return

