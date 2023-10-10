function [x, h] = gmm_sample(gmm, N, seed)
% gmm_sample - sample from a GMM
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2018-09-29
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

% VERSIONS
% 2021-07-21: v0.78 ABC - initial version
%                       - added seed


% set seed for reproducible results
if nargin<3
  seed = 1000;
end
rng(seed, 'twister');

K = length(gmm.pdf);
D = length(gmm.pdf{1}.mean);
  
% sample hidden states
cumprior = cumsum(gmm.prior);
h = multirnd(cumprior, N);

% sample data
x = zeros(N, D);

for j=1:K
  % get the sqrtm of covariance
  tmpC = chol(gmm.pdf{j}.cov, 'lower');
  
  % find points in this cluster
  ii = find(h==j);
  
  % sample
  x(ii,:) = mvgaussrnd(length(ii), gmm.pdf{j}.mean, tmpC); 
end

function x = multirnd(cumprob, n)
% sample a multinomial
f = rand(n, 1);
x = sum(bsxfun(@gt, f, cumprob(:)'),2)+1;

% sample a multivariate Gaussian
function X = mvgaussrnd(N, mn, scov)
D = length(mn);
X = bsxfun(@plus, randn(N, D)*scov', mn(:)');
