function [kld, ll1, ll2, data1] = vbhmm_kld(hmm1, hmm2, data1, opt)
% vbhmm_kld - compute approximation to the KL divergence between two HMMs
%
% KL divergence is a distance measure between probability distributions.
%
%   [kld, ll1, ll2, data1] = vbhmm_kld(hmm1, hmm2, data1, opt)
%
% INPUTS
%   hmm1 = 1st HMM learned with vbhmm_learn
%   hmm2 = 2nd HMM
%  data1 = cell array of fixation sequences (as in vbhmm_learn) used to train hmm1.
%          or [N, T], where N is the number of samples and T the sample-lengths.
%    opt = 'n' - normalize log-likelihood of each sequence by the length of the sequence (default)
%           '' - no special options
%
% OUTPUTS
%   kld = approximate KL divergence between the two HMMs -- a measure of dissimilarity
%         between the two HMMs.  kld=0 when the two HMMs are identical.
%   ll1 = log-likelihood of data1 using HMM1
%   ll2 = log-likelihood of data1 using HMM2
% data1 = data used to compute KL (either the data1 input or randomly sampled data)
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-01-13
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong

% v0.74 - 2018-12-04 - option to randomly sample data, can output the sampled data

if nargin<4
  opt = 'n';
end

% v0.74 - sample data
if ~iscell(data1)
  N = data1(1);
  T = data1(2);
  [~, data1] = vbhmm_random_sample(hmm1, T, N);
end

ll1 = vbhmm_ll(hmm1, data1, opt);
ll2 = vbhmm_ll(hmm2, data1, opt);

kld = mean(ll1-ll2);