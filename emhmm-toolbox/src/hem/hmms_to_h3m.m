function [h3m] = hmms_to_h3m(hmms, covmode)
% hmms_to_h3m - convert a list of HMMs into H3M toolbox format
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
% Copyright (c) 2017-01-13
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong

% VERSIONS
%  2017-01-19: modify for fixation duration
%  2018-11-23: v0.74 - handle empty HMMs w/ component weight=0

% get number of HMMs
h3m.K = length(hmms);

% 2018-11-23: v0.74 - do the normalization later;
%h3m.omega = ones(1,h3m.K) / h3m.K;
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
    s = length(hmms{j}.prior);
    tempHmm.prior = hmms{j}.prior;
    tempHmm.A     = hmms{j}.trans;
    
    % for each emission density, convert to H3M format
    for i = 1:s
      tempHmm.emit{i}.type     = 'gmm';  % a GMM with one component
      tempHmm.emit{i}.nin      = nin;
      tempHmm.emit{i}.ncentres = 1;
      tempHmm.emit{i}.priors   = 1;
      tempHmm.emit{i}.centres  = hmms{j}.pdf{i}.mean;
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


%h3m_old = old_code(hmms);

%% OLD CODE from tim - for checking %%
function largeH3M = old_code(hmms)
largeH3M = {};
for i = 1:10
    largeH3M = ConstructH3M(i,largeH3M,hmms{i,1});
end

function largeH3M = ConstructH3M(id,largeH3M,hmm)
s = size(hmm.prior,1);
for i = 1:s
        tempHmm.emit{1,i}.type = 'gmm';
        tempHmm.emit{1,i}.nin = 2;
        tempHmm.emit{1,i}.ncentres = 1;
        tempHmm.emit{1,i}.covar_type = 'diag';
        tempHmm.emit{1,i}.priors = 1;
        tempHmm.emit{1,i}.centres = hmm.pdf{i,1}.mean;
        tempHmm.emit{1,i}.covars(1,1) = hmm.pdf{i,1}.cov(1,1);
        tempHmm.emit{1,i}.covars(1,2) = hmm.pdf{i,1}.cov(2,2);
end
tempHmm.prior = hmm.prior;
tempHmm.A = hmm.trans;
largeH3M.hmm{1,id} = tempHmm;
largeH3M.K = length(largeH3M.hmm);
omegas = ones(1,largeH3M.K);
largeH3M.omega = omegas;
largeH3M.omega = largeH3M.omega/sum(largeH3M.omega);
