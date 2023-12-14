clear
close all
 
%% options structure for the experiment %%%%%%%%%%%%%%%%%%%%%%%%
OPT = struct;
 
% set DEBUGMODE=1, only the first 10 subjects and first 5 stimuli and used
%                  and only 5 trials are used for VB-HMM and VHEM.
%                  This is mainly to test the code.
% set DEBUGMODE=2, only the first 10 subjects, first 15 stimuli, and use
%                  default number of trials for VB-HMM and VHEM
%                  This is mainly to test the code.
% set DEBUGMODE=0 for the real experiment.
OPT.DEBUGMODE = 0;
 
%% the following data is required:
%
% 1) fixation XLS file - with the following columns:
%      SubjectID - subject name (same as before)
%      TrialID   - trial number (same as before)
%      StimuliID - the name of the stimuli image viewed (required for co-clustering)
%      FixX - fixation X coordinate (same as before)
%      FixY - fixation Y coordinate (same as before)
%
% NOTE: For one subject, TrialID should be unique across all stimuli.
%       I.e., different stimuli cannot have the same TrialID.
%
% 2) imgsize = [xsize, ysize] of stimuli 
%    - we assume that all stimuli images are the same size.
%
% 3) imginfo XLS file - contains the information about each stimuli, 
%    and the valid box region  for that stimuli (e.g., if there is a border around 
%    the actual stimuli). It has the following columns:
%      StimuliID - image name
%      Xlo - lower X coordinates of valid box
%      Xhi - upper X coordinates
%      Ylo - lower Y coordinates
%      Yhi - upper Y coordinates
% 
% 4) imgdir - location of the images
% 5) outdir - where to save the output images and mat files.
 
root = 'H:\OneDrive - The University Of Hong Kong\mscoco';
 
% Notes: the stimuli used in identification and explanation are
% identical. So the stimuli_info files are shared 
 
% Uncleaned EXP Fixation Data
OPT.fixationfile = [root '/cocluster_with_duration_exp_fix_11_18_2023_12_31.xlsx'];
OPT.imgdir       = [root '/Stimuli_Exp/'];
OPT.imgsize      = [1270; 784];
OPT.imginfo      = [root '/image_info_screen.xlsx'];
OPT.outdir       = [root '/exp_results_1group_1203/'];
OPT.HEM_K         = 1;    % fit one group to identify outliers

% individual HMM settings (set vbopt like this)
OPT.S               = 4:10;  % number of ROIs to try   
OPT.vbopt.numtrials = 200;

% co-clustering settings
OPT.vb_cocluster  = 0;
OPT.HEM_S         = [];   % use median
OPT.SWAP_CLUSTERS = []; % swap cluster IDs for consistency w/ paper
OPT.hemopt.seed   = 3000;
OPT.hemopt.trials = 200;
outindmat = 'individual_hmms_1203.mat';
outgrpmat = 'cogroup_hmms_1203.mat';  % file name for group mat

% Edge cases list
edge_cases_path = [root '\edgecase.xlsx'];
edge_cases = readtable(edge_cases_path);

if OPT.vb_cocluster == 0
    OPT.edgecase_crghmm_outdir = [OPT.outdir '/edgecase_hem_crghmm_outdir/']
else
    OPT.edgecase_crghmm_outdir = [OPT.outdir '/edgecase_vbhem_crghmm_outdir/']
end
[success, msg] = mkdir(OPT.edgecase_crghmm_outdir);

do_indhmms = 0;    % run individual HMMs (set to 0 to reuse results)
do_cluster = 1;    % run co-clustering
OPT.DEBUGMODE = 0; % for code testing, set debugging mode to 1
 
%% options for saving images
OPT.SAVE_FIXATION_PLOTS = 0;  % change th
OPT.SAVE_GROUP_MEMBERS_PLOTS = 1; % change this to 1 if you want to plot cluster members with group hmm.
OPT.IMAGE_EXT = 'png';        % format for output images
OPT.SORT_BY_DIFF = 1;         % 1 = when generating image outputs, sort by KL divergence between HMMs
 
OPT
 
%% figure options
figopts = {'Position', [20 0 1240 800]};
stpopts = {0.005, 0.01, 0.01};
subplotwidth = 5;


  
%% load the data and individual HMMs from before %%%%%%%%%%%%%%%%
infile = [OPT.outdir outindmat];
fprintf('loading individual hmms: %s\n', infile);  
if ~exist(infile, 'file')
error('individual hmm mat file does not exist. set do_indhmms=1 and re-run');
end

OPT_new = OPT; % save current OPT
load(infile)
OPT = OPT_new; % set the new OPT

% read in the images, which are not saved in the mat file
StimuliNamesImagesC = cell(length(StimuliNamesC),2);
StimuliNamesImagesC(:,1) = StimuliNamesC;
for j=1:length(StimuliNamesC)
StimuliNamesImagesC{j,2} = imread([OPT.imgdir StimuliNamesC{j}]);
end

% setup output names - use stimuli name or stimuli index
stimname = @(j) ['_' sprintf('%03d', j)];
stimname = @(j) ['_' StimuliNamesC{j}(1:end-4)];
 
 
 
%% exit if we are not doing clustering
if (do_indhmms && ~do_cluster)
  fprintf('-- just computing individual HMMs --\n');
  return
end
 
 
%% co-clustering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (do_cluster)
 
  %% vb co-clustering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if OPT.vb_cocluster
 
    timeVal = tic;
 
    % default parameters
    hemopt = struct;
    vbhemopt.tau   = round(get_median_length(alldataC(:)));
    vbhemopt.m0    = OPT.imgsize/2;
    vbhemopt.W0    = 0.005;
    vbhemopt.initmode  = 'vbh3m'; 
    vbhemopt.learn_hyps = 1;
    vbhemopt.Nv   = 10;
    vbhemopt.seed = 1000;
    
    % add options from OPT structure
    if isfield(OPT, 'vbhemopt')
      fn = fieldnames(OPT.vbhemopt);
      for i = 1:length(fn)
        ff = fn{i};
        vbhemopt.(ff) = OPT.vbhemopt.(ff);
        if isnumeric(vbhemopt.(ff))
          fprintf('-- set vbhemopt.%s = %s\n', ff, mat2str(vbhemopt.(ff)));
        else
          fprintf('-- set vbhemopt.%s = \n', ff)
          vbhemopt.(ff)
        end
      end
    end
    
    if OPT.DEBUGMODE
      fprintf('<strong>*** DEBUG: set trials to 5</strong>\n');
      vbhemopt.trials = 5;
    end
    
    tstart = tic;
    
    telapsed_vbco = toc(tstart)
    
    for i=1:height(edge_cases)
        img_name = [edge_cases{i,'stimuli'}{1} '.png'];
        new_s = edge_cases{i,'roi_number'};
        img_index = find(strcmp(StimuliNamesC,img_name));
        hmm_subject = hmm_subjects(img_index,:);
        
        % run VB co-clustering
        [vbco] = vbhem_cocluster(hmm_subject, OPT.HEM_K, new_s, vbhemopt);
        
        %% save mat file
        outfile = [OPT.edgecase_crghmm_outdir sprintf('cogroup_hmms_%d_%s_HEM_S_%d.mat',img_index, img_name, new_s)];
        fprintf('saving group MAT file: %s\n', outfile);
        save(outfile, 'vbco','vbhemopt','OPT');
 
    end
    
    %% save mat file
    outfile = [OPT.outdir outgrpmat];
    fprintf('saving group MAT file: %s\n', outfile);
    save(outfile, 'vbco','vbhemopt','OPT');
  
    vbhemTime = toc(timeVal);
    
  else
 
    timeVal = tic;
 
    %% coclustering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hemopt.tau = round(get_median_length(alldataC(:)));
    hemopt.seed = 3000;
    hemopt.trials = 50;
    hemopt.initmode  = 'gmmNew';    
        
    % add options from OPT structure
    if isfield(OPT, 'hemopt')
      fn = fieldnames(OPT.hemopt);
      for i = 1:length(fn)
        ff = fn{i};
        hemopt.(ff) = OPT.hemopt.(ff);
        if isnumeric(hemopt.(ff))
          fprintf('-- set hemopt.%s = %s\n', ff, mat2str(hemopt.(ff)));
        else
          fprintf('-- set hemopt.%s = \n', ff)
          hemopt.(ff)
        end
      end
    end
    
    tstart = tic;
    
    %% run co-clustering
    
    for i=1:height(edge_cases)
        img_name = [edge_cases{i,'stimuli'}{1} '.png'];
        new_s = edge_cases{i,'roi_number'};
        img_index = find(strcmp(StimuliNamesC,img_name));
        hmm_subject = hmm_subjects(img_index,:);
        
        % run non-VB co-clustering
        [cogroup_hmms, trials] = vhem_cocluster(hmm_subject, OPT.HEM_K, new_s, hemopt);
        
        %% save mat file
        outfile = [OPT.edgecase_crghmm_outdir sprintf('cogroup_hmms_%d_%s_HEM_S_%d.mat',img_index, img_name, new_s)];
        fprintf('saving group MAT file: %s\n', outfile);
        save(outfile, 'cogroup_hmms','hemopt','OPT'); 
 
    end
    
    
    
    telapsed_co = toc(tstart)
    
%     %% swap clusters (if desired)
%     if ~isempty(OPT.SWAP_CLUSTERS)     
%       cogroup_hmms_old = cogroup_hmms;
%       for i=1:length(cogroup_hmms)
%         cogroup_hmms{i} = vhem_permute_clusters(cogroup_hmms{i}, OPT.SWAP_CLUSTERS);
%       end
% 
%       cogroup_hmms = cogroup_hmms_old
%     end
%     
%     %% save mat file
%     outfile = [OPT.outdir outgrpmat];
%     fprintf('saving group MAT file: %s\n', outfile);
%     save(outfile, 'cogroup_hmms','hemopt','OPT');  
 
    hemTime = toc(timeVal);
  end
     
else
  %% load previous group model
  infile = [OPT.outdir outgrpmat];
  fprintf('loading group file: %s\n', infile);
  
  if ~exist(infile, 'file')
    error('group hmm mat file does not exist. set do_cluster=1 and re-run');
  end
    
  load(infile)
end
 

