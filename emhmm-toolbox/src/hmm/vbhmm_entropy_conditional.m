function H = vbhmm_entropy_conditional(hmm, t)
% vbhmm_entropy_conditional - compute the conditional entropy of an HMM
%
%   H = vbhmm_entropy_conditional(hmm, t)
%
%   Estimate the conditional entropy of an HMM -- the entropy
%   of a fixation at time t, given the fixation at time t-1.
%   i.e., the entropy of p(x_t | x_{t-1}), averaged over p(x_{t-1})
%
%  INPUT
%    hmm = an HMM
%      t = the time-step (fixation number) to use 
%          default is inf, which means use the steady-state distribution
%
%  OUTPUT
%    H = the conditional entropy of the HMM
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2021-07-21
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

% VERSIONS
% 2021-07-21: v0.78 ABC - initial version

if nargin<2
  t = inf;
end
if t<2
  error('t must be >= 2');
end

D = length(hmm.pdf{1}.mean);

NUMSAMPLES = 10000;

if isinf(t)
  % get steady-state joint x_{t-1}, x_t
  %pytt = vbhmm_prob_steadystateobs2(hmm);
  % get steady state marginal
  %pyt = vbhmm_prob_steadystateobs(hmm);
  
  % get steady-state joint of (x_{t-1}, x_t)
  pytt = vbhmm_prob_jointobs2(hmm, inf);
  % get steady-state marginal of x_{t-1}
  pyt = vbhmm_prob_marginalobs(hmm, inf);

else  
  % use time t
  
  % get joint of (x_{t-1}, x_t)
  pytt = vbhmm_prob_jointobs2(hmm, t);
  
  % get marginal of x_{t-1}
  [tmp] = vbhmm_prob_marginalobs(hmm, t-1);
  pyt = tmp{end};
end


% sample
y = gmm_sample(pytt, NUMSAMPLES);

% compute conditional entropy
H = -mean(log(gmm_pdf(pytt, y)) - log(gmm_pdf(pyt, y(:,1:D))));
