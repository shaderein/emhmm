% demo_vb_faces - example of eye gaze analysis for face recognition 
%     uses variational Bayesian HEM (VBHEM) to automatically estimate
%     the number of groups and number of states during clustering.
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2021-07-29
% Antoni B. Chan, Janet H. Hsiao, Lan Hui
% City University of Hong Kong, University of Hong Kong

% 2021-07-29: v0.80 - Hui, ABC - initial version

clear
close all


%% Data files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% names of files used
xlsname = 'demodata.xls';    % Excel File with fixation data
faceimg = 'ave_face120.png';  % average face image

% names of files for saving results
matfile_individual = 'models_demo_vb_faces_individual.mat';
matfile_group      = 'models_demo_vb_faces_group.mat';

% pause after each step in the analysis
%do_pause = 1;
do_pause = 0;  % uncomment this to skip pausing


%% Load data from xls %%%%%%%%%%%%%%%%%%%%%%%%%%
% see the xls file for the format
[data, SubjNames, TrialNames] = read_xls_fixations(xlsname);

% load face image
img0 = imread(faceimg);
imgsize = size(img0);

%% VB Parameters for HMMs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S = 1:3; % automatically select from K=1 to 3
vbopt.alpha0 = 1;
vbopt.mu0    = [imgsize(2); imgsize(1)]/2; 
vbopt.W0     = 0.001;
vbopt.beta0  = 1;
vbopt.v0     = 10;
vbopt.epsilon0 = 1;
vbopt.showplot = 0;     % set this to 1 to see each subject during learning
vbopt.bgimage = faceimg;
vbopt.seed = 100;  % set random state seed for reproducible results.

% Estimate hyperparameters automatically for each individual
% (remove this option to use the above hyps)
vbopt.learn_hyps = 1;  

%% Learn Subjects' HMMs %%%%%%%%%%%%%%%%%%%%%
% estimate HMMs for each individual (same as demo_faces)
[hmms, Ls] = vbhmm_learn_batch(data, S, vbopt);

%% save to individual models to disk
fprintf('saving individual models to %s\n', matfile_individual);
save(matfile_individual, 'hmms', 'faceimg', 'vbopt');
% later on you can use the following command to read back the models:
% load models_demo_vb_faces_individual.mat
 
pause_msg(do_pause);


%% Cluster HMMs using VBHEM %%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('=== Clustering ===\n');

K = 1:5; % automatically select from 1 to 5 clusters
S = 1:3; % automatically select from 1 to 3 ROIs (states)

% initial hyperparameters (similar to vbhmm_learn)
vbhemopt.alpha0   = 1;
vbhemopt.eta0     = 1;
vbhemopt.m0       = vbopt.mu0; 
vbhemopt.W0       = 0.001;
vbhemopt.lambda0  = 1;
vbhemopt.v0       = 10;
vbhemopt.epsilon0 = 1;
vbhemopt.learn_hyps = 1; % automatically estimate the hyperparameters
vbhemopt.seed = 1001; 
vbhemopt.tau = 5;

% perform VB clustering
group_vbhmms = vbhem_cluster(hmms, K, S, vbhemopt); 

%% plot the groups
vhem_plot(group_vbhmms, faceimg);

% plot the groups and cluster members
vhem_plot_clusters(group_vbhmms, hmms, faceimg);

% plot fixations for each group
vhem_plot_fixations(data, group_vbhmms, faceimg);

% plot fixations for each group w/ transition matrix
vhem_plot_fixations(data, group_vbhmms, faceimg, 'c');

% plot the log-likelihoods of each (K,S) pair
figure
plot(K, group_vbhmms.model_LL, 'x-')
tmpleg = cellfun(@(x) sprintf('S=%d',x), num2cell(S), 'uniformoutput', false);
legend(tmpleg) 
xlabel('K');
ylabel('model LL');
grid on
title('LL for different (K,S) pairs');

% show group membership
fprintf('Group membership (indicies): \n');
for j=1:length(group_vbhmms.groups)
  fprintf('  group %d = %s\n', j, mat2str(group_vbhmms.groups{j}));
end

fprintf('Group membership (SubjectID): \n');
for j=1:length(group_vbhmms.groups)
  fprintf('  group %d = ', j)
  fprintf('%s, ', SubjNames{group_vbhmms.groups{j}})
  fprintf('\n');
end

pause_msg(do_pause);


%% Statistical test %%%%%%%%%%%%%%%%%%%%%%%%%%%
% collect data for group 1 and group 2
data1 = data(group_vbhmms.groups{1});
data2 = data(group_vbhmms.groups{2});

fprintf('\n*** t-test ***\n');

% run t-test for hmm1 
[p, info, lld] = stats_ttest(group_vbhmms.hmms{1}, group_vbhmms.hmms{2}, data1);
fprintf('- test group hmm1 different from group hmm2: t(%d)=%0.4g; p=%0.4f; d=%0.3g\n', info.df, info.tstat, p, info.cohen_d);

% run t-test for hmm2
[p, info, lld] = stats_ttest(group_vbhmms.hmms{2}, group_vbhmms.hmms{1}, data2);
fprintf('- test group hmm2 different from group hmm1: t(%d)=%0.4g; p=%0.4f; d=%0.3g\n', info.df, info.tstat, p, info.cohen_d);


%% get mean log-likelihoods (e.g., to use with a correlation test) %%%%%%%%%
fprintf('\n*** get mean LL for each subject ***\n');

[mLL1] = stats_meanll(group_vbhmms.hmms{1}, data);
[mLL2] = stats_meanll(group_vbhmms.hmms{2}, data);

fprintf('- mean LL for each subject under group hmm1:\n');
mLL1
fprintf('- mean LL for each subject under group hmm2:\n');
mLL2

% compute AB scale (e.g., to use with correlation test)
AB = (mLL1 - mLL2) ./ (abs(mLL1) + abs(mLL2));
fprintf('- AB scale values:\n');
AB


%% save group models to file
fprintf('saving group models to %s\n', matfile_group);
save(matfile_group, 'vbhemopt', 'group_vbhmms', 'faceimg');
% later on you can use the following command to read back the models:
% load models_demo_vb_faces_group.mat

