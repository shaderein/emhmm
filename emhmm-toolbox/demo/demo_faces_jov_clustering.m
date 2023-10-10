% demo_faces_jov_clustering - example of HMM clustering
%
% ---
% For each subject, we train an HMM. The subjects' HMMs are clustered
% using VHEM to obtain the common strategies used.
%
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-01-13
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong

% 2017-08-22: updated for v0.7 - auto hyperparameter estimation
%                              - use vbhmm_learn_batch
%                              - save models to mat files
%                              - added sequence length for VHEM
% 2019-02-21: v0.75 - use vbopt.seed to set random state  
%                     NOTE: this is a major change, so you will need to 
%                           update your own scripts that are based on this demo file.

clear
close all


%% Load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load jov_data.mat 
% jov_data contains the data that was used in 
% Chuk T., Chan, A. B., & Hsiao, J. H. (2014). 
% Understanding eye movements in face recognition using
% hidden Markov models. Journal of Vision 14(8).
% doi:10.1167/14.11.8.

% data is stored in a cell array
% data{i}         = i-th subject
% data{i}{j}      = ... j-th trial
% data{i}{j}(t,:) = ... [x y] location of t-th fixation 

% the number of subjects
N = length(data);

faceimg = 'ave_face120.png';


%% VB Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K = 1:3;  % automatically select from K=1 to 3
vbopt.alpha0 = 1;
vbopt.mu0    = [160;210];  
vbopt.W0     = 0.001;
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


% names of files for saving results
matfile_individual = 'models_demo_faces_jov_clustering_individual.mat';
matfile_group      = 'models_demo_faces_jov_clustering_group.mat';

%% Learn Subject's HMMs %%%%%%%%%%%%%%%%%%%%%

% estimate an HMM for each individual
[hmms, Ls] = vbhmm_learn_batch(data, K, vbopt);

% plot all subjects
for i=1:N
  if mod(i,16)==1
    figure
  end
  subplot(4,4,mod(i-1,16)+1)
  vbhmm_plot_compact(hmms{i}, faceimg);
  title(sprintf('Subject %d', i));
end

%% save to individual models to disk
fprintf('saving individual models to %s\n', matfile_individual);
save(matfile_individual, 'hmms', 'faceimg', 'vbopt');
% later on you can use the following command to read back the models:
% load models_demo_faces_jov_clustering_individual.mat


%% Run HEM clustering (1 cluster) %%%%%%%%%%%%%%%%
fprintf('=== Clustering 1 ===\n');
hemopt.tau = get_median_length(data);  % set the virtual sequence length to the median data length.
hemopt.seed = 1001;  % set random seed

[group_hmms1] = vhem_cluster(hmms, 1, 3, hemopt)

% plot the groups
vhem_plot(group_hmms1, faceimg);


%% Run HEM Clustering (2 clusters) %%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('=== Clustering 2 ===\n');
[group_hmms2] = vhem_cluster(hmms, 2, 3, hemopt)

% plot the groups
vhem_plot(group_hmms2, faceimg);

% plot the groups and cluster members
vhem_plot_clusters(group_hmms2, hmms, faceimg);

% plot fixations for each group
vhem_plot_fixations(data, group_hmms2, faceimg);

% plot fixations for each group w/ transition matrix
vhem_plot_fixations(data, group_hmms2, faceimg, 'c');

% show group membership
fprintf('Group membership: \n');
for j=1:length(group_hmms2.groups)
  fprintf('  group %d = %s\n', j, mat2str(group_hmms2.groups{j}));
end


%% Run HEM Clustering (3 clusters) %%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('=== Clustering 3 ===\n');
[group_hmms3] = vhem_cluster(hmms, 3, 3, hemopt)

% plot the groups
vhem_plot(group_hmms3, faceimg);

% plot the groups and cluster members
vhem_plot_clusters(group_hmms3, hmms, faceimg);

% plot fixations for each group
vhem_plot_fixations(data, group_hmms3, faceimg)

% plot fixations for each group w/ transition matrix
vhem_plot_fixations(data, group_hmms3, faceimg, 'c')

% show group membership
fprintf('Group membership: \n');
for j=1:length(group_hmms3.groups)
  fprintf('  group %d = %s\n', j, mat2str(group_hmms3.groups{j}));
end


% save group models to file
fprintf('saving group models to %s\n', matfile_group);
save(matfile_group, 'hemopt', 'group_hmms1', 'group_hmms2', 'group_hmms3', 'faceimg');
% later on you can use the following command to read back the models:
% load models_demo_faces_jov_clustering_group.mat


