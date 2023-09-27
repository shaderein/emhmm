% demo_conversion - a demo showing how to convert data collected from one face image 
%  to another face image with different size. Then calculate LL 
%  with respect to a model trained on the other face images.
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-05-25
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong

% v0.75 - use vbopt.seed for random state
% v0.78 - fix bug: use old data
clear
%close all


%% Load our data %%%%%%%%%%%%%%%%%%%%%%%%%%
load demodata_old.mat 

% the data is stored in a cell array:
% data{i}         = i-th subject
% data{i}{j}      = ... j-th trial
% data{i}{j}(t,:) = ... [x y] location of t-th fixation 

% the number of subjects
N = length(data);

% our data's face image
our_faceimg = 'face.jpg';


%% Load the model from JOV %%%%%%%%%%%%%%%%%%%%%%
% (these were saved from demo_face_jov_clustering.m)
model = load('jov_models.mat');

% show the model
vhem_plot(model.group_hmms2, model.faceimg);


%% setup the coordinate conversion %%%%%%%%%%%%%%
% we normalize based on the midpoint between the eyes, and the distance from
% this midpoint to the mouth level.

% for our images:
our_C  = [258, 186];  % midpoint between eyes
our_D  = 121;         % distance from eyes midpoint to mouth level

% for the other images from JOV:
other_C = [160,195];  % midpoint between eyes
other_D = 90;         % distance from eyes midpoint to mouth level

% show a plot
plot_face_coord('new data', our_faceimg, our_C, our_D, ...
                'JOV model', model.faceimg, other_C, other_D);

%% do the coordinate conversion %%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_s = convert_face_coord(data, our_C, our_D, other_C, other_D);

% show the original data and converted data
figure
for i=1:min(3,N) % just show 3
  subplot(2,3,i)
  plot_fixations(data{i}, our_faceimg, [], '')
  title(sprintf('Subject %d (original)', i));
  subplot(2,3,3+i)
  plot_fixations(data_s{i}, model.faceimg, [], '');
  title(sprintf('Subject %d (converted)', i));
end


%% calculate Log-likelihood of converted data with the JOV model %%
fprintf('\n*** get mean LL for each subject ***\n');

[mLL1] = stats_meanll(model.group_hmms2.hmms{1}, data_s);
[mLL2] = stats_meanll(model.group_hmms2.hmms{2}, data_s);

fprintf('- mean LL for each subject under hmm1:\n');
mLL1
fprintf('- mean LL for each subject under hmm2:\n');
mLL2

%% assign to a group HMM %%%%
fprintf('\n*** subject assignments to HMMs ***\n');
[~, label] = max([mLL1;mLL2], [], 1);
label

fprintf('\n*** groups ***\n');
for i=1:2
  groups{i} = find(label==i);
  fprintf('group %d: ', i);
  fprintf('%d ', groups{i});
  fprintf('\n');
end
