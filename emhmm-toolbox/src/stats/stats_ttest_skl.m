function [p, info, lld] = stats_ttest_skl(hmms1, hmms2, data1, data2, opt)
% stats_ttest_skl - run t-test to compare two HMMs using symmetric KL divergence
%
%  [p, info] = stats_ttest_skl(hmms1, hmms2, data1, data2, opt)
%
%   The LL difference of hmm1(data1)-hmm2(data1) and
%                        hmm2(data2)-hmm1(data2)
%   are concatenated together. The average of this vector is the symmetric
%   KL divergence between hmm1 and hmm2.  The t-test indicates the 
%   symmetric KL > 0.
%
% INPUTS
%   hmms1 = 1st HMM (or a cell array of HMMS, one for each subject)
%   hmms2 = 2nd HMM (or a cell array of HMMS, one for each subject)
%   data1 = cell array of data associated with each subject in hmms1
%           data1{i}{j} - for each subject i (belonging to hmms1{i}), and j-th trial
%   data2 = cell array of data associated with each subject in hmms2
%    opt = 'n' - normalize log-likelihood of each sequence by the length of the sequence (default)
%           '' - no special options
%
% OUTPUTS
%      p = the p-value
%   info = other info from the t-test (Cohen's d, df, etc)
%    lld = the SKL values used for the t-test
%
% SEE ALSO
%   vbhmm_kld
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2021-04-04
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

% 2021-04-04: ABC - initial version v0.78

if nargin<5
  opt = 'n';
end

% compute LL of HMM1(data1) - HMM2(data1)
[~,~,lld12] = stats_ttest(hmms1, hmms2, data1, opt);
% compute LL of HMM2(data2) - HMM1(data2)
[~,~,lld21] = stats_ttest(hmms2, hmms1, data2, opt);

% concatenate -- taking the average will give the SKL
lld = [lld12,lld21];

% do a t-test: since average(lld) >= 0.  
% so we can use a right-tailed test to check if the mean is > 0.
[h,p,ci,stats] = ttest(lld, 0, 'tail', 'right');

% compute Cohen's d for one-sample t-test
d = mean(lld) / std(lld);

info = stats;
info.ci = ci;
info.skl_avg  = 0.5*(mean(lld12)+mean(lld21));
info.skl      = mean(lld);
info.cohen_d  = d;

return

% test code
%   ttKL = []
%   ttp = []
%   for j=1:Nstimuli
%     [p,info,lld] = stats_ttest_skl(cogroup_hmms{j}.hmms{1}, cogroup_hmms{j}.hmms{2}, alldataC(j,grps{1}), alldataC(j,grps{2}));
%     ttKL(j) = info.skl;
%     ttp(j) = p;
%   end
%   figure(1)
%   plot(KLsum, 'b')
%   hold on
%   plot(ttKL, 'r')
%   hold off
%  
