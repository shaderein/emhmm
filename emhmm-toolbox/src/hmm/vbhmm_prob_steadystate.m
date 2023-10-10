function p = vbhmm_prob_steadystate(hmm)
% vbhmm_prob_steadystate - compute the steady-state probability of the ROIs (hidden states)
%
% p = vbhmm_prob_steadystate(hmm)
%
%  INPUT:  hmm - HMM
%  OUTPUT: p   - probability vector
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2018-06-13
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

% v0.73 - initial version

if isfield(hmm, 'group_map')
  for j=1:length(hmm.trans)
    p{j} = computep(hmm.trans{j}, hmm.prior{j});
  end
else
  p = computep(hmm.trans, hmm.prior);
end

% do the calculation
% https://nicolewhite.github.io/2014/06/10/steady-state-transition-matrix.html
%
% need to solve p = Ap
function p = computep(A, prior)
d = size(A,1);

% check if transition matrix is degenerate (no transitions)
if all(all(abs(A-eye(d)) < 1e-6))
  % use the prior as steady-state
  p = prior;
else
  p = [(A'-eye(d,d)); ones(1,d)] \ [zeros(d,1);1];
end
