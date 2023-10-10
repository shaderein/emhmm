% demo_faces_jov_compare - example of comparing HMMs
%
% ---
% For each subject, we separated trials that led to correct responses
% and trials that led to incorrect (wrong) responses.
% A separate HMM is learned for correct trials and wrong trials.
% Then, we compare the "correct" HMM and "wrong" HMM.
%
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-01-13
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong

% 2017-08-22: updated for v0.7 - auto hyperparameter estimation
%                              - use vbhmm_learn_batch
%                              - save models to mat files
% 2019-02-21: v0.75 - use vbopt.seed to set random state  
%                     NOTE: this is a major change, so you will need to 
%                           update your own scripts that are based on this demo file.

clear
close all


%% Load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load jov_data.mat 
% jov_data contains the data that was used for 
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
K = 1:3; 
vbopt.alpha0 = 1;
vbopt.mu0    = [160;210]; 
vbopt.W0     = 0.001; 
vbopt.beta0  = 1; 
vbopt.v0     = 10; 
vbopt.epsilon0 = 1;
vbopt.showplot = 0;
vbopt.bgimage = faceimg;

vbopt.seed = 1000;  % set random state for reproducible results.

% names for saving results
matfile_individual = 'models_demo_jov_compare_individual.mat';


%% Learn Subject's HMMs %%%%%%%%%%%%%%%%%%%%%

% train one subject at a time
for i=1:N
  % share hyperparameters between correct and wrong trials of the same subject
  % (assumes that the Subject's ROIs are similar between correct/wrong trials)
  vbopt.learn_hyps_batch = 1;
  
  % setup the batch containing correct/wrong trials for one subject
  data_tmp = {data_correct{i}, data_wrong{i}};  
  [hmms_tmp, Ls] = vbhmm_learn_batch(data_tmp, K, vbopt);
  
  % save
  hmms_correct{i} = hmms_tmp{1};
  hmms_wrong{i}   = hmms_tmp{2};  
  
  % view correct/wrong models
  figure(100)
  clf
  subplot(1,2,1)
  vbhmm_plot_compact(hmms_correct{i}, faceimg);
  title(sprintf('Subject %d - correct', i));
  subplot(1,2,2)
  vbhmm_plot_compact(hmms_wrong{i}, faceimg);
  title(sprintf('Subject %d - wrong', i));
  drawnow
end

if 0
  % this is another way to do it...
  %
  % assume that correct & wrong trials of the same subject should have 
  % different hyperparameters. I.e., the Subject's ROIs are different between
  % correct/wrong trials.
  vbopt.learn_hyps = 1;

  % collect correct and wrong trials into one batch
  data_all = {data_correct{:}, data_wrong{:}};
  
  % learn an HMMs for each subject and correct/wrong trials
  [hmms_all, Ls] = vbhmm_learn_batch(data_all, K, vbopt);
  
  % separate into correct and wrong HMMs
  hmms_correct = hmms_all(1:N);
  hmms_wrong = hmms_all(N+(1:N));
  
  %% view correct/wrong models
  for i=1:N
    figure(100)
    clf
    subplot(1,2,1)
    vbhmm_plot_compact(hmms_correct{i}, faceimg);
    title(sprintf('Subject %d - correct', i));
    subplot(1,2,2)
    vbhmm_plot_compact(hmms_wrong{i}, faceimg);
    title(sprintf('Subject %d - wrong', i));
    
    pause
  end
end


%% save to individual models to disk
fprintf('saving individual models to %s\n', matfile_individual);
save(matfile_individual, 'hmms_correct', 'hmms_wrong', 'faceimg', 'vbopt');
% later on you can use the following command to read back the models:
% load models_demo_faces_jov_compare_individual.mat


%% Run statistical tests %%
% see if correct HMMs are different from wrong HMMs.
fprintf('=== correct vs. wrong ===\n');
[p, info, lld] = stats_ttest(hmms_correct, hmms_wrong, data_correct);
p
info

fprintf('=== wrong vs. correct ===\n');
[p, info, lld] = stats_ttest(hmms_wrong, hmms_correct, data_wrong);
p 
info


return

%% compare with old
for i=1:N
  figure(100)
  clf
  subplot(2,2,1)
  vbhmm_plot_compact(hmms_correct{i}, faceimg);
  title(sprintf('Subject %d - correct (LL=%0.4g)', i, hmms_correct{i}.LL));
  subplot(2,2,2)
  vbhmm_plot_compact(hmms_wrong{i}, faceimg);
  title(sprintf('Subject %d - wrong (LL=%0.4g)', i, hmms_wrong{i}.LL));
  
  subplot(2,2,3)
  vbhmm_plot_compact(OLD.hmms_correct{i}, faceimg);
  title(sprintf('OLD Subject %d - correct (LL=%0.4g)', i, OLD.hmms_correct{i}.LL));
  subplot(2,2,4)
  vbhmm_plot_compact(OLD.hmms_wrong{i}, faceimg);
  title(sprintf('OLD Subject %d - wrong (LL=%0.4g)', i, OLD.hmms_wrong{i}.LL));
  pause
end


