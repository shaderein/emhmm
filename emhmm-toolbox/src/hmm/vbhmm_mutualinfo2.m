function MI = vbhmm_mutualinfo2(hmm, t)
% vbhmm_mutualinfo2 - compute the mutual information between x_{t-1} and x_t in the HMM
%
%   MI = vbhmm_mutualinfo2(hmm, t)
%
%   Estimate the mutual information of fixations t and t-1 of an HMM.
%
%  INPUT
%    hmm = an HMM
%      t = the time-step (fixation number). 
%
%  OUTPUT
%    MI = the mutual information
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2018-10-25
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

% VERSIONS
% v0.74 - initial version

if t<2
  error('t must be >= 2');
end

NUMSAMPLES = 10000;
D = length(hmm.pdf{1}.mean);

% use time t

% get marginal and entropy
pxm = vbhmm_prob_marginalobs(hmm, t);
pyt  = pxm{t};
pyt1 = pxm{t-1};

% get joint of (x_{t-1}, x_t)
pytt = vbhmm_prob_jointobs2(hmm, t);

% sample
y = gmm_sample(pytt, NUMSAMPLES);

% compute mutual information
MI = mean( log(gmm_pdf(pytt, y)) ...
         - log(gmm_pdf(pyt1, y(:,1:D))) ...
         - log(gmm_pdf(pyt,  y(:,(D+1):end))) ...
         );

