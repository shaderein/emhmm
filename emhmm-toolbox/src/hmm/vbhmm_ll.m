function [loglik, errors] = vbhmm_ll(hmm, data, opt)
% vbhmm_ll - compute the log-likelihood of sequences under an HMM
%
%  [loglik, errors] = vbhmm_ll(hmm, data, opt)
%
% INPUTS
%   hmm = HMM learned with vbhmm_learn
%  data = cell array of fixation sequences (as in data passed to vbhmm_learn)
%   opt = 'n' - normalize log-likelihood of each sequence by the length of the sequence (default) --
%               longer sequences typically have lower log-likelihoods (log-likelihood for each
%               fixation is negative, and they add). Normalization will make log-likelihoods
%               of different-length sequences more comparable.
%          '' - no normalization
%
% OUTPUTS
%   loglik = log-likelihoods of each sequence [1 x N]
%   errors = whether errors occurred for each sequence [1 x N]
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-01-13
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong

% 2017-05-17: BUG FIX - Sequence normalization was wrong for length-1 sequences.
%                       The old code erroneously thought these were length-2 sequences.
%                       This should only affect comparison of log-likelihoods (for correlation tests).
%

% NOTE: will not work with groups

if nargin<3
  opt = 'n';
end

if ~iscell(data)
  data = {data};
end

prior = hmm.prior;
transmat = hmm.trans;

K = length(hmm.pdf);
N = length(data);

fb_qnorm = zeros(1,N);

% for each sequence
for n = 1:N  
  tdata = data{n};
  tT = size(tdata,1);
  
  % calculate observation likelihoods
  logrho = zeros(tT,K);  
  for k=1:K
    logrho(:,k) = log(mvnpdf(tdata, hmm.pdf{k}.mean, hmm.pdf{k}.cov));
  end
  
  % forward algorithm
  t_logPiTilde = prior(:)';
  t_logATilde = transmat;
  t_logrho = exp(logrho);
  
  t_alpha = []; t_beta = []; t_c = [];
  
  if tT >= 1
    %forward
    t_alpha(1,:) = t_logPiTilde.*t_logrho(1,:);
    
    % 2016-04-29 ABC: rescale for numerical stability (otherwise values get too small)
    t_c(1) = sum(t_alpha(1,:));
    t_alpha(1,:) = t_alpha(1,:) / t_c(1);
    
    if tT > 1
      for i=2:tT
        t_alpha(i,:) = t_alpha(i-1,:)*t_logATilde.*t_logrho(i,:);
        
        % 2016-04-29 ABC: rescale for numerical stability
        t_c(i) = sum(t_alpha(i,:));
        t_alpha(i,:) = t_alpha(i,:) / t_c(i);
      end;
    end
    
  end
  
  % from scaling constants
  fb_qnorm(n) = sum(log(t_c));  
end

loglik = fb_qnorm;

% normalize
if any(opt=='n')
  %loglik = loglik ./ reshape(cellfun('length', data), size(loglik));
   
  % 2017-05-17: BUG fix -- make the sequence length correct when only one fixation in a trial
  [tmpnr, tmpnc] = cellfun(@size, data);
  loglik = loglik ./ reshape(tmpnr, size(loglik));
end

errors = isinf(loglik);
