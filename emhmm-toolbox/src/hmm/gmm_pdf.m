function [P] = gmm_pdf(gmm, X)
% gmm_pdf - compute the likelihood under a GMM
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2018-09-29
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

% VERSIONS
% 2021-07-21: v0.78 - ABC - initial version

% compute probabilities
P = zeros(size(X,1),1);

for k=1:length(gmm.pdf)
    myp = mvnpdf(X, gmm.pdf{k}.mean, gmm.pdf{k}.cov);
    P = P + gmm.prior(k)*myp;
end
