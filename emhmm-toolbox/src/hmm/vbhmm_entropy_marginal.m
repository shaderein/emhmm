function H = vbhmm_entropy_marginal(hmm, T)
% vbhmm_entropy_marginal - compute the entropy of marginal distributions
%
%   H = vbhmm_entropy_marginal(hmm, T)
%
%  Compute the entropy of the marginal distributions of an HMM, 
%  i.e.  the entropy of fixations at time 1, time 2, etc.
%  
% INPUTS
%   hmm = an HMM
%     T = time length (how many time steps to compute)
%
% OUTPUTS
%   H = [1 x T] the entropy of the marginal distribution p(x_t) for time 1:T
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2021-07-21
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

% VERSIONS
% 2021-07-21: v0.78 ABC - initial version

% get marginals
pxm = vbhmm_prob_marginalobs(hmm, T);

% compute entropy via sampling
NUMSAMPLES = 10000;
for i=1:T
  H(i) = gmm_entropy(pxm{i}, NUMSAMPLES);
end


