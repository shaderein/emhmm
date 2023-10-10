% demo_vb_faces_jov_duration - example of HMMs with fixation location and duration
%  using VBHEM for clustering, which automatically selects the number of 
% groups and number of states.
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2021-07-29
% Antoni B. Chan, Janet H. Hsiao, Lan Hui
% City University of Hong Kong, University of Hong Kong

% 2021-07-29: v0.80 - Hui - initial version


clear
close all

%% Load data with fixation location and duration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[data,SubjectIDs,TrialIDs] = read_xls_fixations('jov_duration.xls');

% data is stored in a cell array
% data{i}         = i-th subject
% data{i}{j}      = ... j-th trial
% data{i}{j}(t,:) = ... [x y d] location of t-th fixation. duration "d" is in milliseconds 

% data is also in this mat file:
%load jov_duration.mat 


% the number of subjects
%N = length(data);

% only look at 10 in this demo
N = 15;  
data = data(1:10);

faceimg = 'ave_face120.png';

% to be faster for live demo
%vbopt.numtrials = 5;  
%vbhemopt.trials = 5; 

%% VB Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K = 2:3;
vbopt.alpha0 = 1;
vbopt.mu0    = [160;210;250];  % the image center, and 250ms duration
vbopt.W0     = [0.001, 0.001, 0.0016]; % stdev of 31 pixels for ROIs, and stddev of 25ms for duration
vbopt.beta0  = 1;
vbopt.v0     = 10; 
vbopt.epsilon0 = 1;
vbopt.showplot = 1;     % show each subject during learning
vbopt.bgimage = faceimg;

vbopt.seed = 1000;  % set random state for reproducible results.

% Estimate hyperparameters automatically for each individual
% (remove this option to use the above hyps)
vbopt.learn_hyps = 1;  

% Estimate hyperparameters automatically as a batch (all individuals share the same hyperparameters)
% (uncomment below to enable the option)
%vbopt.learn_hyps_batch = 1;


%% Learn Subject's HMMs %%%%%%%%%%%%%%%%%%%%%

% estimate an HMM for each individual
[hmms, Ls] = vbhmm_learn_batch(data, K, vbopt);


% plot one subject
vbhmm_plot(hmms{6}, data{6}, faceimg);
figure; vbhmm_plot_compact(hmms{6}, faceimg);


%% Run VBHEM clustering %%%%%%%%%%%%%%%%

K = 1:5; % select among 1-5 groups
S = 1:4; % select among 1-4 ROIs

% initial hyperparameters (similar to vbhmm_learn)
vbhemopt.alpha0 = 10;
vbhemopt.eta0   = 1;
vbhemopt.m0     = [160;210;250]; 
vbhemopt.W0     = [0.001, 0.001, 0.0016]; 
vbhemopt.lambda0  = 1;   % is beta0 in vbhmm_learn
vbhemopt.v0       = 5;
vbhemopt.epsilon0 = 1;
vbhemopt.learn_hyps = 1; % automatically estimate the hyperparameters
vbhemopt.seed = 1800;  
vbhemopt.tau= 5;
vbhemopt.trials = 50;

% run vbhem clustring
[group_vbhmms] = vbhem_cluster(hmms, K, S, vbhemopt);

%% plot the models
vhem_plot(group_vbhmms, faceimg);

% plot fixations for each group
vhem_plot_fixations(data, group_vbhmms, faceimg);

% plot fixations for each group w/ transition
vhem_plot_fixations(data, group_vbhmms, faceimg, 'c');

% show group membership
fprintf('Group membership: \n');
for j=1:length(group_vbhmms.groups)
  fprintf('  group %d = %s\n', j, mat2str(group_vbhmms.groups{j}));
end

% save group models to file
% names of files for saving results
matfile_group      = 'models_demo_vb_faces_duration_group.mat';

fprintf('saving group models to %s\n', matfile_group);
save(matfile_group, 'vbhemopt',  'group_vbhmms', 'faceimg');
% later on you can use the following command to read back the models:
% load models_demo_faces_duration_group.mat


