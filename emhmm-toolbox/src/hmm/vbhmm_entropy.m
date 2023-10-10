function [ent, nll] = vbhmm_entropy(hmm, data)
% vbhmm_entropy - compute approximation to the entropy of an HMM given data
%
% Entropy is a measure of certainty/uncertainty. Higher values means more uncertain or random, 
% while Entropy of 0 means completely deterministic.
%
%   [ent, nll] = vbhmm_entropy(hmm, data)
%
% INPUTS
%   hmm = HMM learned with vbhmm_learn
%  data = cell array of fixation sequences (as in vbhmm_learn) used to train hmm
%
% OUTPUTS
%   ent = the approximate entropy of the HMM. 
%   nll  = negative log-likelihood of data using HMM.
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2018-09-05
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

% v0.74 - initial version

% compute negative likelihoods, and normalize by sequence length
nll = -vbhmm_ll(hmm, data, 'n');

% compute entropy
ent = mean(nll);