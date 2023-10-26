% Author: Alice Yang.
% Date: 12/12/2022
% Email: yumeng@g.ucla.edu
% HKU_ORIB_Eye_Movement_Post_Model_Data_Processing
%%% this script will analyze the AB scale and also perform KL divergence
%%% analysis.....

%%%%%
%%% Input: 1) vbcogroup_hmms_veh_0814.mat from EMHMM
%%%%%%%%%% 2) cluster_asd_useful_variables_analysis_veh_hmm_08_15_2023_09_00.mat
%%%%%%%%%% 3) data_ready_2analyze_acc_char_entropy_rt_08_22_2023_11_48.csv
%%% Output: 1) 
% =========================================================================
clc
clear

%% setting of the script
taskType = 2; % CHANGE 1 -> vehicle detection task.
saveKLD = 0; % 1 -> save, 0 -> no save
subgroup = "no"; % if you only want to extract particular ppts
processed_data_name = "data_ready_2analyze_acc_char_entropy_rt_08_22_2023_11_48.csv";
%% set up path and load data
% set up the path for the code to run, make sure the fixation xlsx is in
% the same folder as this script
myname = mfilename('fullpath'); % the current m file's fullpath
[filepath name ext] = fileparts(myname); % get the current path and file name

% go to one level above the current folder
idcs = strfind(filepath,'/'); % find the "/", so that you can go to one level above
parent_dir = filepath(1:idcs(end)-1); % only include path before the last "/"
addpath(genpath(parent_dir))

% create a save folder for the results
save_dir = fullfile(parent_dir,"/processed_data/"); % where to save
disp("This is the saving directory: " + save_dir)

% load the files to get the groups

% import the files you need
if taskType == 1
    % load('vbcogroup_hmms_veh_0814.mat')
    load('cluster_asd_useful_variables_diffAB_model_veh_hmm_09_03_2023_17_30.mat')
    % load('cluster_asd_useful_variables_analysis_veh_hmm_08_15_2023_09_00.mat')
    img_info = readtable("autism_object_detection_emhmm/veh_id_stimuli_info.xlsx"); % load the image's information
    sub_sheet = readtable("clean_fix_data/with_duration_detection_veh_fix_08_14_2023_22_35.xlsx"); % load the fixation sheet.
else
    % load('vbcogroup_hmms_hum_0814.mat')
    load('cluster_asd_useful_variables_diffAB_model_hum_hmm_09_03_2023_17_48.mat')
    % load('cluster_asd_useful_variables_analysis_hum_hmm_08_15_2023_18_24.mat')
    img_info = readtable("autism_object_detection_emhmm/hum_id_stimuli_info.xlsx");
    sub_sheet = readtable("clean_fix_data/with_duration_detection_hum_fix_08_14_2023_23_22.xlsx");
end

% if you want to save a particular group of participants 
if subgroup == "yes"
    sub_info = unique(sub_sheet.SubjectID);
    sub_info = sub_info(1:49);
else
    sub_info = unique(sub_sheet.SubjectID);
end

%% find the group info and put the sub info with the group and AB scale
% group1 = vbco.cogroup_hmms{1,1}.groups{1,1}'; % grab the group info from struct
% group2 = vbco.cogroup_hmms{1,1}.groups{1,2}';
% group = zeros(numel(sub_info),1); % pre allocate a cell
% group(group1',:) = 1;
% group(group2',:) = 2;
grp1 = find(AB>0); 
grp2 = find(AB<0);
group = zeros(numel(sub_info),1); % pre allocate a cell
group(grp1) = 1;
group(grp2) = 2;

group_info = table(sub_info,group,AB);
group_info_ab = table(sub_info,group,AB);


for i=1:size(group_info_ab,1)
    IsASD = contains(group_info_ab.sub_info{i},"A");
    % if subjectId has A, then the partiStatus is ASD
    if IsASD == 1
        group_info_ab.partiStatus{i} = 'ASD';
    else
        group_info_ab.partiStatus{i} = 'Control';
    end
end

group_info_ab1 = extractBeforeAccordingAnotherCol(group_info_ab,1,"Veh",'subjectNo'); % extract the subject number for each participant
group_info_ab1.subjectNo = str2double(group_info_ab1.subjectNo);
% import the processed data 
processed_data = readtable(processed_data_name);
processed_data = renamevars(processed_data,["subjectNo_behavioralData","partiStatus_behavioralData"],["subjectNo","partiStatus"]);
data_with_ab = outerjoin(group_info_ab1,processed_data,"Keys",["subjectNo","partiStatus"]);
data_with_ab2 = removevars(data_with_ab,["subjectNo_processed_data","subjectNo_entropyData","partiStatus_processed_data","partiStatus_entropyData"]);
data_with_ab2 = renamevars(data_with_ab2,["partiStatus_group_info_ab1","subjectNo_group_info_ab1"],["partiStatus","subjectNo"]);

%% save the grp info and ab results in csv file
% crate the naming for the files and save 
csvname = sprintf('5data_wh_ab_%s.csv',datestr(now,'mm_dd_yyyy_HH_MM'));
writetable(data_with_ab2,fullfile(save_dir,csvname));
disp('saving the data into csv files here: ' + save_dir)

%% create the KLD for ANOVA analysis
do_KLD = input("Do you want to conduct KLD analysis to examine whether the two EM group differed from each other? \n")
if do_KLD ==1 
LL1_grp1_with_grp = [LL1_grp1 ones(size(LL1_grp1,1),1)];    % subj in group 1 under group 1 model
LL2_grp1_with_grp = [LL2_grp1 ones(size(LL2_grp2,1),1).*2]; % subj in group 2 under group 1 model
LL1_grp2_with_grp = [LL1_grp2 ones(size(LL1_grp1,1),1)];    % subj in group 1 under group 2 model
LL2_grp2_with_grp = [LL2_grp2 ones(size(LL2_grp2,1),1).*2]; % subj in group 2 under group 2 model
model1 = [LL1_grp1_with_grp;LL2_grp1_with_grp];
model2 = [LL1_grp2_with_grp; LL2_grp2_with_grp];
all_subj_2model = [model1 model2];
all_subj_2model = array2table(all_subj_2model,'VariableNames', ...
    {'LL_model1','EM_group','LL_model2','EM_group_2'}); % create the table
subj_2model_t = array2table([all_subj_2model.EM_group ...
    all_subj_2model.LL_model1 all_subj_2model.LL_model2], ...
    'VariableNames',{'EM_group','LL_model1','LL_model2'}); % re-format the table

%%%%%%%%%%% performance repeated mea ANOVA, model x group anova %%%%%%%%%%%
em_model_level = table([1 2]','VariableNames',{'em_model'});
rm_em_models = fitrm(subj_2model_t,'LL_model1-LL_model2~EM_group','WithinDesign',em_model_level);
disp('Please see the repeated measure anova table below for modelx group...')
rep_mea_anova_4_model = ranova(rm_em_models) % conduct repeated mea anova for model x group

%%%%%%% calculate each stimulus's average LL, stimulus x group anova %%%%%%
LL1_original = LL(:,:,1); % get LL from model 1 for each subject for 160 images
LL2_original = LL(:,:,2);
mean_LL_by_img = array2table((LL1_original + LL2_original)./2); % calculate
% the average LL for each subj for 160 images

group_t = array2table(group);  % convert the group membership to a table
mean_LL_t = [group_t mean_LL_by_img]; % mean LL and group membership table
image_level = array2table((1:160)',"VariableNames",{'stimulus'});
rm = fitrm(mean_LL_t,'Var1-Var160~group','WithinDesign',image_level);
disp('Please see the repeated measure anova table below for stimulus x group...')
rep_mea_anova_4_stimulus = ranova(rm) % this is the repeated measure anova

mean_LL_group_name = mean_LL_t.Properties.VariableNames(:,2:161);
mean_LL_within_mat = table2array(mean_LL_t(:,2:161));
[tbl,rm] = simple_mixed_anova(mean_LL_within_mat,mean_LL_t.group, {'stimulus'}, {'group'})

else
    disp("I did not do KLD analysis. Please input 1 if you want me to. Thanks! ")
end

% %% save the data for KLD for both models
% if saveKLD == 1
%     disp('The script is saving KLD estimate. ')
%     if taskType == 1
%         writetable(subj_2model_t,fullfile(filepath,"veh_KLD_estimate_sheet_4_anova_300_rerun.csv"));
%     else
%         writetable(subj_2model_t,fullfile(filepath,"hum_KLD_estimate_sheet_4_anova.csv"));
%     end
% 
% else
%     disp('The script did not save KLD estimate. ')
% end


