%% HKU_ORIB_Coclustering_Fixation_Data
% =========================================================================
% Goal: run the raw excel sheets from dataviewer, usually named as 
% "Veh_IdTask_60Parti_Fixation1207.xlsx"
% The output file is also an excel sheet, the format is as below
% SubjectID, TrialID, StimuliID, FixX, FixY
% version: Alice Yang     '06-Dec-2022 20:01:23'
% =========================================================================

clc
clear

% set up the path for the code to run, make sure the fixation xlsx is in
% the same folder as this script
myname = mfilename('fullpath'); % the current m file's fullpath
[filepath name ext] = fileparts(myname); % get the path and current file's name

% create the folder directory and cd it as the current folder
addpath(genpath(filepath)); % add all the folders under this directory


%% import the raw files 
% input the fixations from the raw file 
fixInputR = readtable("../bdd/fixation/raw_fixation_report_jinhan_231002/human_exp_raw.xlsx"); % CHANGE the file's name
fixInput1 = array2table([fixInputR.CURRENT_FIX_X fixInputR.CURRENT_FIX_Y ...
    fixInputR.CURRENT_FIX_DURATION fixInputR.TRIAL_INDEX]);
fixInput1.Properties.VariableNames = {'CURRENT_FIX_X','CURRENT_FIX_Y','CURRENT_FIX_DURATION','Trial'}; % rename the variables
% fixInput1.CURRENT_FIX_Y = fixInput1.CURRENT_FIX_Y - 30; % normalize the y coordinates so that the fixations will match

% for the pilot, the session label is number, so i add this part to combine
% all the information 
fixInput2 = array2table([fixInputR.image fixInputR.RECORDING_SESSION_LABEL]);
fixInput2.Properties.VariableNames = {'Image','SessionLabel'};

% uncomment the script below! if the session label becomes a string!!!
fixInput2 = array2table([fixInputR.image fixInputR.RECORDING_SESSION_LABEL]);
fixInput2.Properties.VariableNames = {'Image','SessionLabel'};

fixInput = [fixInput1 fixInput2]; % concat the tables into one

%% create the fixation file to run coclustering 
cocluster_data = removevars(fixInput,"CURRENT_FIX_DURATION");
check   = contains(cocluster_data.Image, '_id.jpg'); % check to see if the rows has contained practice trials
cocluster_data = cocluster_data(~check, :          ); 

% Count trials in each block
% human
block_start = [0,49,98,146]; % [49,49,48,14];

% vehicle
block_start = [0,49,98,148]; % [49,49,50,12];


% Explanation data grouped into blocks in the format 001ETx
% remove block suffix and 
% use cumulative trial numbering

% % Extract the relevant columns from the table
sessionLabels = cocluster_data.SessionLabel;
trials = cocluster_data.Trial;

% Get the last character of each session label
blocks = cellfun(@(x) str2double(x(end)), sessionLabels);

% Remove the last character from each session label
SubjectID = cellfun(@(x) x(1:end-1), sessionLabels, 'UniformOutput', false);

% Vectorize the block_start function
blockStartValues = arrayfun(@(x) block_start(x), blocks);

% Add the corresponding block_start value to each trial
trials = trials + blockStartValues;

TrialID = num2str(trials);
num = size(TrialID);
trial = repmat('Trial  ',num(1),1);
TrialID = strcat(trial,TrialID);
StimuliID = cocluster_data.Image; 
FixX = cocluster_data.CURRENT_FIX_X;
FixY = cocluster_data.CURRENT_FIX_Y;
cocluster_data = table(SubjectID, TrialID, StimuliID, FixX, FixY);

% save the export file 
filename = sprintf('../bdd/fixation/hum_exp_screen_%s.xlsx', datestr(now,'mm_dd_yyyy_HH_MM'));
writetable(cocluster_data, [filepath '/' filename]);