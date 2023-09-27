function H = gmm_entropy(gmm, NumSamples)
% gmm_entropy - compute the entropy of GMM
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2018-09-29
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

% VERSIONS
% 2021-07-21: v0.78 - ABC - initial version


if nargin<2
  NumSamples = 10000;
end

X = gmm_sample(gmm, NumSamples);

H = -mean(log(gmm_pdf(gmm, X)));