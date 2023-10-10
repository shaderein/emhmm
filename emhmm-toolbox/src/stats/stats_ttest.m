function [p, info, lld] = stats_ttest(hmms1, hmms2, data1, opt)
% stats_ttest - run t-test to compare two HMMs
%
%  [p, info, lld] = stats_ttest(hmms1, hmms2, data1, opt)
%
% INPUTS
%   hmms1 = 1st HMM (or a cell array of HMMS, one for each subject)
%   hmms2 = 2nd HMM (or a cell array of HMMS, one for each subject)
%   data1 = cell array of data associated with each subject in hmms1
%           data1{i}{j} - for each subject i (belonging to hmms1{i}), and j-th trial
%    opt = 'n' - normalize log-likelihood of each sequence by the length of the sequence (default)
%           '' - no special options
%
% OUTPUTS
%      p = the p-value
%   info = other info from the t-test (df, t-statistic, confidence interval, Cohen's d, etc)
%    lld = the log-likelihood differences used for the t-test
%
% SEE ALSO
%   vbhmm_kld
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-01-13
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong

% 2017-05-25: ABC - update to use stats_meanll
% 2021-04-09: ABC - skip empty data
% 2021-07-21: ABC - calculate Cohen's d, returned in info

if nargin<4
  opt = 'n';
end

N = length(data1);
allLL1 = [];
allLL2 = [];

% for each subject
for i=1:N
  % get the HMMs
  if iscell(hmms1)
    myhmm1 = hmms1{i};
  else
    myhmm1 = hmms1;
  end
  if iscell(hmms2)
    myhmm2 = hmms2{i};
  else
    myhmm2 = hmms2;
  end
  
  % 2021-04-09: ABC - skip empty data
  if ~isempty(data1{i})    
    % calculate log-likelihoods under each HMM
    ll1 = vbhmm_ll(myhmm1, data1{i}, opt);
    ll2 = vbhmm_ll(myhmm2, data1{i}, opt);
    
    % calculate the mean log-likelihood for this subject
    allLL1(end+1) = mean(ll1);
    allLL2(end+1) = mean(ll2);
  end
end

% do a t-test between the mean log-likelihoods of all subjects
% theoretically, ll1 >= ll2, since hmm1 was learned from data1.
% in other words, (ll1-ll2) >= 0.  
% so we can use a right-tailed test to check if the mean is > 0.
[h,p,ci,stats] = ttest(allLL1, allLL2, 'tail', 'right');

% tim's version
% [h,p,ci,stats] = ttest(allLL1, allLL2, 0.05, 'right');

% v0.78 - compute Cohen's d
%d = computeCohen_d(allLL1, allLL2, 'paired');
d = mean(allLL1-allLL2) / std(allLL1-allLL2);


% output results
lld = allLL1 - allLL2;
info = stats;
info.ci = ci;
info.cohen_d = d;


