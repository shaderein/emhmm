function [Agroup_hmms] = inivb_vhemco(hmms_subjects,vbhemopt,K,S)
% initial vbhem co-clustering via vhem co-clustering
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2021-07-29
% Antoni B. Chan, Janet H. Hsiao, Lan Hui
% City University of Hong Kong, University of Hong Kong

% 2021-07-29: v0.80 - Hui - initial version

sizeall = size(hmms_subjects);
stimulisize = sizeall(1);
subjectsize = sizeall(2);

hemopt.seed = vbhemopt.seed;
hemopt.trials = vbhemopt.trials;
hemopt.keepTopK = vbhemopt.trials; % must less or equal to hemopt.trials
hemopt.tau = vbhemopt.tau;
hemopt.initmode ='gmmNew';
hemopt.sortclusters =vbhemopt.sortclusters;

[~,h3m_out_trials] = vhem_cocluster(hmms_subjects, K, S, hemopt);

Agroup_hmms = cell(hemopt.keepTopK,1);
for j =1:hemopt.keepTopK
    h3m_out2 = h3m_out_trials{1, 1}.topK_h3ms{1, j};
    grouphmms = [];
    
    for i=1:stimulisize
        %make it compatible with LogL and LogLs
        h3m_out2.h3m_rs(i).LogLs = h3m_out2.h3m_rs(i).LogL;
        grouphmm = h3m_to_hmms(h3m_out2.h3m_rs(i));
        grouphmms = [grouphmms; grouphmm];
    end
    

    % sort clusters
    group_hmms = cell(stimulisize,1);
    if ~isempty(hemopt.sortclusters)
        for i=1:stimulisize
            group_hmms{i} = vbhmm_standardize(grouphmms(i), hemopt.sortclusters);
            group_hmms{i}.hemopt = hemopt;
        end
    end
    Agroup_hmms{j} = group_hmms;
    
end

