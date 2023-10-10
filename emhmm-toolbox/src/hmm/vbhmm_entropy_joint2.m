function [ent] = vbhmm_entropy_joint2(hmm, t)
% vbhmm_entropy_joint - compute joint entropy of HMM for sequence pairs (x_{t-1}, x_t)
%
%   [ent, nll] = vbhmm_entropy_joint2(hmm, t)
%
%  computes entropy of fixations at time t-1 and time t from the HMM,
%  i.e. entropy of (x_{t-1},x_t)
%
% INPUTS
%   hmm = HMM learned with vbhmm_learn
%     t = time step (fixtion number)
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

if t<2
  error('t must be >= 2');
end

% get joint
pxj = vbhmm_prob_jointobs2(hmm, t);

% compute entropy via sampling
NUMSAMPLES = 10000;
ent = gmm_entropy(pxj, NUMSAMPLES);


