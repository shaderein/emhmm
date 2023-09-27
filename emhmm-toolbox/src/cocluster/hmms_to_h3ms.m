function [h3m] = hmms_to_h3ms(hmms, covmode)
% hmms_to_h3m - convert a list of HMMs into H3M toolbox format (co-clustering)
%
%     h3m = hmms_to_h3m(hmms, covmode)
%
% INPUT:  hmms    = two dimentional cell array of HMMs {stimuli}{subject}
%         covmode = 'diag' - use diagonal covariance
%                 = 'full' - use full covariance
% OUTPUT: h3ms  = HMM mixture format for H3M toolbox (for co-clustering)
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2019-07-13
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

% VERSIONS
%  2017-01-19: modify for fixation duration
%  2020-03-12: v0.77 - initial version


% get number of HMMs
sizeall = size(hmms)
h3m.Kstimuli = sizeall(1);
h3m.Ksubject = sizeall(2);
%h3m.omega = ones(1,h3m.K) / h3m.K;

% input dimension
nin = length(hmms{1,1}.pdf{1}.mean);

% convert each HMM to the H3M toolbox format
for i=1:h3m.Kstimuli
    for j=1:h3m.Ksubject 
      clear tempHmm

      % state prior, transition matrix 
      tempHmm.prior = hmms{i,j}.prior;
      tempHmm.A     = hmms{i,j}.trans;
      s = length(tempHmm.prior);

      % for each emission density, convert to H3M format
      for k = 1:s
        tempHmm.emit{k}.type     = 'gmm';  % a GMM with one component
        tempHmm.emit{k}.nin      = nin;
        tempHmm.emit{k}.ncentres = 1;
        tempHmm.emit{k}.priors   = 1;
        tempHmm.emit{k}.centres  = hmms{i,j}.pdf{k}.mean;
        switch(covmode)
          case 'diag'
            tempHmm.emit{k}.covar_type = 'diag';        
            tempHmm.emit{k}.covars = diag(hmms{i,j}.pdf{k}.cov)';
          case 'full'
            tempHmm.emit{k}.covar_type = 'full';
            tempHmm.emit{k}.covars     = hmms{i,j}.pdf{k}.cov;
          otherwise
            error('unknown covmode')
        end
      end

      % save to new structure
      h3m.hmm{i,j} = tempHmm;
    end
end

return