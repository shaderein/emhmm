function [h, data] = vbhmm_random_sample(hmm, T, N)
% vbhmm_random_sample - generate random samples from an HMM
%
%   [h, data] = vbhmm_random_sample(hmm, T, N)
%
%  INPUTS
%   hmm = input HMM
%     T = [scalar value] - the length of each sequence
%         [1xL] vector ] - sample a length from this distribution, 
%                          where T(j) = probability of length=j
%     N = number of sequences to sample
%
%  OUTPUT
%    h = 1xN cell array of hidden states
% data = 1xN cell array of fixation data
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2018-09-06
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

% 2018-09-06 - v0.74 initial version

h = cell(1,N);

cumprior = cumsum(hmm.prior);
cumtrans = cumsum(hmm.trans, 2);

if ~isscalar(T)
  cumT = cumsum(T);
end

for n=1:N
  % get sequence length
  if isscalar(T)
    myT = T;
  else
    % sample T
    myT = multirnd(cumT);
  end
  
  % sample hidden states
  hh = zeros(myT,1);
  hh(1) = multirnd(cumprior);
  for t=2:myT
    hh(t) = multirnd(cumtrans(hh(t-1),:));
  end
  h{n} = hh;
end

% sample data
if (nargout > 1)
  D = length(hmm.pdf{1}.mean);
  K = length(hmm.pdf);
  
  data = cell(1,N);
  for n=1:N
    myT = length(h{n});
    data{n} = zeros(myT,D);
  end
  
  for j=1:K
    % get the sqrtm of covariance
    tmpC = chol(hmm.pdf{j}.cov, 'lower');
    
    % sample
    for n=1:N
      ii = find(h{n}==j);
      data{n}(ii,:) = mvgaussrnd(length(ii), hmm.pdf{j}.mean, tmpC);
    end
  end  
end


function x = multirnd(cumprob)
% sample a multinomial
f = rand(1);
x = sum(f>cumprob)+1;

% sample a multivariate Gaussian
function X = mvgaussrnd(N, mn, scov)
D = length(mn);
X = bsxfun(@plus, randn(N, D)*scov', mn(:)');
