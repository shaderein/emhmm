function [p] = vbhmm_prob_marginalstate(hmm, T)
% vbhmm_prob_marginalstate - compute marginal probability of hidden states
%
%   [p] = vbhmm_prob_marginalstate(hmm, T)
%
%   INPUT:  hmm - the hmm
%           T   - time length
%
%   OUTPUT: p = [S x T] - the t-th column is the marginal probability
%                         distribution for the state at time t.
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2021-07-21
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

% VERSIONS
% 2021-07-21: v0.78 - ABC - initial version


S = length(hmm.pdf);
p = zeros(S, T);

p(:, 1) = hmm.prior;

for t=2:T
   p(:, t) =  hmm.trans'*p(:, t-1);
end

