function [pr] = vbhmm_prob_state(hmm, stateseq)
% vbhmm_prob_state - calculate the probability of a hidden-state sequence
%
% [pr] = vbhmm_prob_state(hmm, stateseq)
%
% INPUTS
%       hmm = HMM learned with vbhmm_learn
%  stateseq = a single state sequence, [1 x L] vector;
%             or a cell array of state sequences: stateseq{i} is the i-th state sequence
%
% OUTPUTS
%   pr = probability of the state sequence;
%        or probability of the state sequences [1 x N]
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-05-25
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong

% v0.73 - rename to vbhmm_prob_state

if ~iscell(stateseq)
  stateseq = {stateseq};
end

prior = hmm.prior;
transmat = hmm.trans;

K = length(hmm.pdf);
N = length(stateseq);

pr = zeros(1,N);

% for each seq
for i=1:N
  X = stateseq{i};
  
  % prior
  mypr = prior(X(1));
  
  % transitions
  for t=2:length(X)
    mypr = mypr * transmat(X(t-1), X(t));
  end
  
  pr(i) = mypr;
end

