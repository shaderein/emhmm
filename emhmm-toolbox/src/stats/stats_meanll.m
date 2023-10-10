function [LLs] = stats_meanll(hmm, data, opt)
% stats_meanll - calculate the mean log-likelihood of subjects' data under an hmm
%
%  [LLs] = stats_meanll(hmms, data, opt)
%
% INPUTS
%   hmm  = the HMM
%   data = cell array of multiple subjects' data
%          data1{i}{j} - for each subject i, and j-th trial
%    opt = 'n' - normalize log-likelihood of each sequence by the length of the sequence (default)
%           '' - no special options
%
% OUTPUTS
%    LLs = vector of mean log-likelihoods [1 x N]
%
% SEE ALSO
%   stats_ttest
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-01-13
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong

if nargin<3
  opt = 'n';
end

N = length(data);
LLs = zeros(1,N);

% for each subject
for i=1:N  
  % calculate log-likelihoods under each HMM
  ll1 = vbhmm_ll(hmm, data{i}, opt);
  
  % calculate the mean log-likelihood for this subject
  LLs(i) = mean(ll1);
end
