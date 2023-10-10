function [ent] = vbhmm_entropy_joint(hmm, T)
% vbhmm_entropy_joint - compute joint entropy of HMM for sequence 1:T
%
%   [ent, nll] = vbhmm_entropy_joint(hmm, T)
%
%  computes entropy of length-T sequences from the HMM, i.e., x_1, ..., x_T
%
% INPUTS
%   hmm = HMM learned with vbhmm_learn
%     T = sequence length
%
% OUTPUTS
%   ent = the approximate entropy of the HMM. 
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2018-10-25
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

% VERSIONS
% v0.74 - initial version

NUMSAMPLES = 1000;

[h, data] = vbhmm_random_sample(hmm, T, NUMSAMPLES);

% compute negative likelihoods, and normalize by sequence length
nll = -vbhmm_ll(hmm, data, '');

% compute entropy
ent = mean(nll);
