% demo_PBR_model - demo showing how to use the Holistic/Analytic models from the PBR paper.
%
% The PBR paper is:
% Cynthia Y.H. Chan, Antoni B. Chan, Tatia M.C. Lee, and Janet H. Hsiao,
% "Eye Movement Patterns in Face Recognition are Associated with Cognitive Decline
% in Older Adults.", Psychonomic Bulletin & Review, to appear 2018.
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2018-05-15
% Antoni B. Chan, Janet H. Hsiao, Cynthia Chan
% City University of Hong Kong, University of Hong Kong

% VERSIONS
%   2018-05-16: v0.72 - initial version
%   2019-02-21: v0.75 - update for vbopt.seed
clear
close all



%% Load our data %%%%%%%%%%%%%%%%%%%%%%%%%%

% select which dataset to load
datamode = 1

switch(datamode)
  case 1
    % some example data
    load Example_HMM_PBR.mat
    data = Example_HMM_PBR(:,3);
    
    % our data's face image
    our_faceimg = '../demo/ave_face120.png';
    
    % we normalize based on the midpoint between the eyes, and the distance from
    % this midpoint to the mouth level.
    our_C = [160,195];  % midpoint between eyes (EDIT IF NECESSARY)
    our_D = 90;         % distance from eyes midpoint to mouth level (EDIT IF NECESSARY)
    
  case 2
    % all JOV data
    load ../demo/jov_data.mat
    
    % our data's face image
    our_faceimg = '../demo/ave_face120.png';
    
    % we normalize based on the midpoint between the eyes, and the distance from
    % this midpoint to the mouth level.
    our_C = [160,195];  % midpoint between eyes (EDIT IF NECESSARY)
    our_D = 90;         % distance from eyes midpoint to mouth level (EDIT IF NECESSARY)
    
  case 3
    % DEMO data
    
    %load demodata.mat
    [data, SubjNames, TrialNames] = read_xls_fixations('../demo/demodata_old.xls');
    
    % our data's face image
    our_faceimg = '../demo/face.jpg';
    
    % we normalize based on the midpoint between the eyes, and the distance from
    % this midpoint to the mouth level.
    our_C  = [258, 186];  % midpoint between eyes (EDIT IF NECESSARY)
    our_D  = 121;         % distance from eyes midpoint to mouth level (EDIT IF NECESSARY)
end

% the number of subjects
N = length(data);

%% Load the representative models %%%%%%%%%%%%%%%%%%%%%%

% Load the models from PBR
model = load('PBR_models.mat');
modelname = 'PBR';

% show the model, please check if they are the models you expected
vhem_plot(model.group_hmms2, model.faceimg, 'c', [], model.modelName);


%% do the coordinate conversion %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Please check whether the midpoint between eyes and the eye-to-mouth distance are
% correct
plot_face_coord('our data', our_faceimg, our_C, our_D, ...
                'PBR data', model.faceimg, model.other_C, model.other_D);

% do the actual conversion
data_s = convert_face_coord(data, our_C, our_D, model.other_C, model.other_D);

% show the original data and converted data for the first 3 participants
% Please check whether the  fixations are converted properly
figure
for i=1:min(3,N) % just show 3
  subplot(2,3,i)
  plot_fixations(data{i}, our_faceimg, [], '')
  title(sprintf('Subject %d (original)', i));
  subplot(2,3,3+i)
  plot_fixations(data_s{i}, model.faceimg, [], '');
  title(sprintf('Subject %d (converted)', i));
end

%% calculate Log-likelihood of converted data with the model %%
% The output will be stored in LL_all under workspace

fprintf('\n*** get mean LL for each subject ***\n');

LL_all = [];
for i=1:length(model.group_hmms2.hmms)
    [mLL] = stats_meanll(model.group_hmms2.hmms{i}, data_s);
    LL_all = [LL_all;mLL];
    fprintf('Column %d refers to the LL of the %s pattern\n',i, model.modelName{i});
end
LL_all=LL_all'


%% compute H-A scale
fprintf('\n*** H-A scale ***\n');
HAscale = (LL_all(:,1) - LL_all(:,2)) ./ sum(abs(LL_all), 2)

% sort HA values in ascending order
[HAscale_sort, inds] = sort(HAscale);

% select 8 values
inds2 = inds(round(linspace(1,length(inds),8)));

%% make a nice plot
figure
% plot the models
vhem_plot(model.group_hmms2, model.faceimg, 'c', [2,3,3,1], model.modelName);
% plot the scale
subplot(2,3,2)
plot(HAscale, zeros(size(HAscale)), '.');
hold on
plot(HAscale(inds2), zeros(size(inds2)), 'o');
hold off
grid on
xlabel(sprintf('<--- %s -----              ----- %s --->', model.modelName{2}, model.modelName{1}));
title('HA scale');
xmm = max(abs(HAscale));
axis([-xmm, xmm, -1, 1]);

% make a plot of the 8 subjects
for i=1:length(inds2)
  subplot(2,8,8+i)
  ii = inds2(i);
  plot_fixations(data_s{ii}, model.faceimg, [], '')
  title(sprintf('HA = %0.4f', HAscale(ii)));
end


%% assign to a group HMM %%%%
fprintf('\n*** subject assignments to HMMs ***\n');
[~, label] = max(LL_all', [], 1);

fprintf('\n*** groups ***\n');
for i=1:size(model.group_hmms2.hmms,2)
  groups{i} = find(label==i);
  fprintf('group %d - %s: ', i, model.modelName{i});
  fprintf('%d ', groups{i});
  fprintf('\n');
end


