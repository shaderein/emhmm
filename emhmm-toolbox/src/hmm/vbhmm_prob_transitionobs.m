function [pdfx] = vbhmm_prob_transitionobs(hmm, xt1, t)
% vbhmm_prob_transitionobs - compute the conditional transition probability of an HMM
%
%   [pdf] = vbhmm_prob_transitionobs(hmm, xt1, t)
%
%  Compute the transition probability p(x_t | x_{t-1})
%
%   INPUT:  hmm - the hmm
%           xt1 - [P x D] conditioned observation x_{t-1}. Each row is one observation (D-dim vector). 
%                 there are P observations
%             t - the time step (fixation number)
%                 [default = inf, which means steady-state]
%
%   OUTPUT: pdf = {1 x P} - cell array of GMMs, one for each conditioned x_{t-1}
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2018-09-29
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

% VERSIONS
% 2021-07-21: v0.78 - ABC - initial version

if nargin<3
  t = inf;
end
if t<2
  error('t must be >= 2');
end


P = size(xt1, 1);
S = length(hmm.pdf);

if isinf(t)
  % get the steady state probabilites
  pz = vbhmm_prob_steadystate(hmm);
else
  % get state marginal at t-1
  tmp = vbhmm_prob_marginalstate(hmm, t-1);
  pz = tmp(:,end);
end


% posterior p(z_{t-1} | x_{t-1})

% initialize with px
pzt1xt1 = repmat(pz(:), [1 P]);

% likelihood of previous xt1 | zt1
for k=1:S
    myp = mvnpdf(xt1, hmm.pdf{k}.mean, hmm.pdf{k}.cov);
    pzt1xt1(k,:) = pzt1xt1(k,:) * myp;
end

% normalize
pzt1xt1 = bsxfun(@times, pzt1xt1, 1./sum(pzt1xt1,1));

% p(z_t | x_{t-1})
pztxt1 = hmm.trans'*pzt1xt1;

% now make the final probability as a Gaussian mixture
pdfx = cell(1,P);

for i=1:P
   % construct GMM for each t
   pdfx{i}.prior = pztxt1(:,i);
   pdfx{i}.pdf   = hmm.pdf;
end


