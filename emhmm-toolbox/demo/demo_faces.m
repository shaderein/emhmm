% demo_faces - example of eye gaze analysis for face recognition 
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2019-06-14
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong

% 2017-08-22: updated for v0.7 - auto hyperparameter estimation
%                              - use vbhmm_learn_batch
%                              - save models to mat files
%                              - added sequence length for VHEM
% 2019-01-12: v0.75 - use vbopt.seed and hemopt.seed to set random state for reproducible results
%                     NOTE: this is a major change, so you will need to 
%                           update your own scripts that are based on this demo file.
% 2019-06-14: v0.76 - clean-up script for CogSci2019 tutorial, new dataset
% 2021-07-21: v0.78 - print out Cohen's d
%                   - use median number of states for HEM

clear
close all

%% Data files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% names of files used
xlsname = 'demodata.xls';    % Excel File with fixation data
faceimg = 'ave_face120.png';  % average face image

% names of files for saving results
matfile_individual = 'models_demo_faces_individual.mat';
matfile_group      = 'models_demo_faces_group.mat';


% pause after each step in the analysis
do_pause = 1;   
%do_pause = 0;  % uncomment this to skip pausing


%% Load data from xls %%%%%%%%%%%%%%%%%%%%%%%%%%
% see the xls file for the format
[data, SubjNames, TrialNames] = read_xls_fixations(xlsname);

% the data is read and separated by subject and trial, and stored in a cell array:
% data{i}         = i-th subject
% data{i}{j}      = ... j-th trial
% data{i}{j}(t,:) = ... [x y] location of t-th fixation 

% the same data is stored in a mat file.
%load demodata.mat 

% the number of subjects
N = length(data);

% load image
img0 = imread(faceimg);
imgsize = size(img0);

%% VB Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K = 1:3;      % automatically select from K=1 to 3
vbopt.alpha0 = 1;
vbopt.mu0    = [imgsize(2); imgsize(1)]/2;  % center of the image
vbopt.W0     = 0.001;
vbopt.beta0  = 1;
vbopt.v0     = 10; 
vbopt.epsilon0 = 1;
vbopt.showplot = 0;     % show each subject during learning
vbopt.bgimage = faceimg;

vbopt.seed = 1000;  % set random state seed for reproducible results.

% Estimate hyperparameters automatically for each individual.
% the values specified above are used as the initial values.
% (remove this option to use the above hyps)
vbopt.learn_hyps = 1;  

% Estimate hyperparameters automatically as a batch (all individuals share the same hyperparameters)
% (uncomment below to enable this option)
%vbopt.learn_hyps_batch = 1;

% Normally the ROIs are learned from the data. 
% To specify a set of fixed ROIs, uncomment the two lines below. 
% See the documentation FAQ for more details.
%vbopt.fixed_rois = {[120,195,80,40], [202,195,80,40], ...
%                    [159,237,40,80], [159,285,80,40]}; 
%K = length(vbopt.fixed_rois);


%% Learn Subjects' HMMs %%%%%%%%%%%%%%%%%%%%%

% estimate HMMs for each individual
[hmms, Ls] = vbhmm_learn_batch(data, K, vbopt);

% show subject 4
mys = 4;
vbhmm_plot(hmms{mys}, data{mys}, faceimg);

% compact plot
figure, vbhmm_plot_compact(hmms{mys}, faceimg);

% show fixations and compact plot
figure
subplot(2,1,1)
plot_fixations(data{mys}, faceimg, [], 's');
subplot(2,1,2)
vbhmm_plot_compact(hmms{mys}, faceimg, 'r', data{mys});

pause_msg(do_pause);

% plot each subject
figure
for i=1:N
  subtightplot(2,5,i)
  vbhmm_plot_compact(hmms{i}, faceimg);
  title(sprintf('SubjectID=%s (index=%d)', SubjNames{i}, i));
end

%% save to individual models to disk
fprintf('saving individual models to %s\n', matfile_individual);
save(matfile_individual, 'hmms', 'faceimg', 'SubjNames', 'vbopt');
% later on you can use the following command to read back the models:
% load models_demo_faces_individual.mat

pause_msg(do_pause);


%% Run HEM clustering (1 group) %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% summarize all subjects with one HMM
fprintf('=== Clustering (1 group) ===\n');
hemopt.tau = get_median_length(data);  % set the virtual sequence length to the median data length.
hemopt.seed = 1001;  % set random state seed for reproducible results.

% number of states in HMMs for HEM
% [] means use the median number of states in the input HMMs
myS = [];

[all_hmms1] = vhem_cluster(hmms, 1, myS, hemopt)  % 1 group

% plot the overall HMM
vhem_plot(all_hmms1, faceimg);

pause_msg(do_pause);


%% Run HEM Clustering (2 groups) %%%%%%%%%%%%%%%%%%%%%%%%%%
% cluster subjects into 2 groups
fprintf('=== Clustering (2 groups) ===\n');
[group_hmms2] = vhem_cluster(hmms, 2, myS, hemopt)  % 2 groups

% plot the groups
vhem_plot(group_hmms2, faceimg);

% plot the groups and cluster members
vhem_plot_clusters(group_hmms2, hmms, faceimg);

% plot fixations for each group
vhem_plot_fixations(data, group_hmms2, faceimg);

% plot fixations for each group w/ transition matrix
vhem_plot_fixations(data, group_hmms2, faceimg, 'c');




% show group membership
fprintf('Group membership (indicies): \n');
for j=1:length(group_hmms2.groups)
  fprintf('  group %d = %s\n', j, mat2str(group_hmms2.groups{j}));
end

fprintf('Group membership (SubjectID): \n');
for j=1:length(group_hmms2.groups)
  fprintf('  group %d = ', j)
  fprintf('%s, ', SubjNames{group_hmms2.groups{j}})
  fprintf('\n');
end

pause_msg(do_pause);

%% get the most likely ROI sequences %%%%%%%%%
fprintf('\n*** top-5 most probable ROI sequences ***\n');
[seqs{1}, seqs_pr{1}] = stats_stateseq(group_hmms2.hmms{1}, 3);
[seqs{2}, seqs_pr{2}] = stats_stateseq(group_hmms2.hmms{2}, 3);

% show the sequences
fprintf('prob   : hmm1 ROI seq \n');
fprintf('-----------------\n');
for i=1:5
  fprintf('%0.4f : ', seqs_pr{1}(i));
  fprintf('%d ',     seqs{1}{i});
  fprintf('\n');
end

% show the sequences
fprintf('\nprob   : hmm2 ROI seq \n');
fprintf('-----------------\n');
for i=1:5
  fprintf('%0.4f : ', seqs_pr{2}(i));
  fprintf('%d ',     seqs{2}{i});
  fprintf('\n');
end

pause_msg(do_pause);


%% Statistical test %%%%%%%%%%%%%%%%%%%%%%%%%%%
% collect data for group 1 and group 2
data1 = data(group_hmms2.groups{1});
data2 = data(group_hmms2.groups{2});

fprintf('\n*** t-test ***\n');

% run t-test for hmm1 
[p, info, lld] = stats_ttest(group_hmms2.hmms{1}, group_hmms2.hmms{2}, data1);
fprintf('- test group hmm1 different from group hmm2: t(%d)=%0.4g; p=%0.4f; d=%0.3g\n', info.df, info.tstat, p, info.cohen_d);

% run t-test for hmm2
[p, info, lld] = stats_ttest(group_hmms2.hmms{2}, group_hmms2.hmms{1}, data2);
fprintf('- test group hmm2 different from group hmm1: t(%d)=%0.4g; p=%0.4f; d=%0.3g\n', info.df, info.tstat, p, info.cohen_d);


%% get mean log-likelihoods (e.g., to use with a correlation test) %%%%%%%%%
fprintf('\n*** get mean LL for each subject ***\n');

[mLL1] = stats_meanll(group_hmms2.hmms{1}, data);
[mLL2] = stats_meanll(group_hmms2.hmms{2}, data);

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
save(matfile_group, 'hemopt', 'all_hmms1', 'group_hmms2', 'faceimg');
% later on you can use the following command to read back the models:
% load models_demo_faces_group.mat


