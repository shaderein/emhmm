function [pdfx] = vbhmm_prob_marginalobs(hmm, T)
% vbhmm_prob_marginalobs - compute marginal probability of observations of an HMM
%
%   [pdf] = vbhmm_prob_marginalobs(hmm, T)
%
%  Computes the marginal distribution p(x_t), for t=1:T
%
%   INPUT:  hmm - the hmm
%           T   - time length
%                 inf = just return steady-state marginal distribution
%
%   OUTPUT: pdf = {1 x T} - cell array of GMMs, one for each p(x_t)
%               = GMM, if T=inf
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2018-09-29
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

% VERSIONS
% 2021-07-21: v0.78 - ABC - initial version

if isinf(T)
  % steady-state
  
  % get the steady-state state probabilites
  pz = vbhmm_prob_steadystate(hmm);  
  
  pdfx.prior = pz(:);
  pdfx.pdf   = hmm.pdf;
else
  % get the marginal state probabilites
  pz = vbhmm_prob_marginalstate(hmm, T);
    
  pdfx = cell(1,T);
  for t=1:T
    % construct GMM for each t
    pdfx{t}.prior = pz(:,t);
    pdfx{t}.pdf   = hmm.pdf;
  end
end
