# EMHMM Analysis

## `./bdd/`

Data and Results folder for EMHMM runs

## `./helpers/`

Helper functions. Remember to add this folder to your path.

## `./scripts/`

1. `fixations_preprocessing.m` adapted from Alice's `Convert_Fixation_2_Cocluster_Format.m` with snippets to map data in different blocks to the single corresponding participant.

2. `cocluster_categorization.m` the main co-clustering script adapted from `demo_cocluster_Categorization_final.m` from emhmm-toolbox and Jennifer's script. Perform 1-group/2-group co-clustering using different cases.

3. `remove_outliers.m` adapted from `vhem_plot_fixations.m` in the emhmm-toolbox that can be called within a 1-group co-clustering run using `cocluster_categorization.m`

4. `replace_group_hmms.m`. This script should be run after obtaining the cogroup hmms for specific edge case stimuli. Set the file name and index of the edge case among all stimuli. It will replace the corresponding parameters in the cogroup hmms of all cases.

### `./scripts/Jennier/`

Original co-clustering script used by Jennifer. Semantic analysis example from Jennifer.

### `./scripts/Alice`

Scripts shared by Alice that are not yet used.

# XAI

1. Script of generate human attention maps for similarity analysis: `H:\Projects\HKU_XAI_Project\XAI_Similarity_1\Main_gen_human_saliency_map_2.m`

# Paths

```
% Human
   
%Categorization
% OPT.fixationfile = 'bdd/fixation/hum_id_fix_12_17_2022_17_58.xlsx';
% OPT.imgdir       = 'bdd/images/orib_hum_id_task_resized/';
% OPT.imgsize      = [1024; 576];
% OPT.imginfo      = 'bdd/image_info/hum_id_stimuli_info.xlsx';
% OPT.outdir       = 'bdd/results/identification/ORIB-data-vb-results-lab-hum-local/';

%Explanation
% OPT.fixationfile = 'bdd/fixation/hum_exp_fix_10_04_2023_02_54.xlsx';
% OPT.imgdir       = 'bdd/images/orib_hum_id_task_resized/';
% OPT.imgsize      = [1024; 576];
% OPT.imginfo      = 'bdd/image_info/hum_id_stimuli_info.xlsx';
% OPT.outdir       = 'bdd/results/explanation/human/';

% Explanation whole screen
% OPT.fixationfile = 'bdd/fixation/hum_exp_screen_10_07_2023_23_49.xlsx';
% OPT.imgdir       = 'bdd/images/human_with_screen/';
% OPT.imgsize      = [1024; 768];
% OPT.imginfo      = 'bdd/image_info/hum_exp_screen_stimuli_info.xlsx';
% OPT.outdir       = 'bdd/results/explanation/231015_human_whole_screen_1-group/';
% OPT.HEM_K         = 1;    % fit one group to identify outliers

% Vehicle

%Categorization
% OPT.fixationfile = 'bdd/fixation/veh_id_fix_12_07_2022_11_00.xlsx';
% OPT.imgdir       = 'bdd/images/orib_veh_id_task_resized/';
% OPT.imgsize      = [1024; 576];
% OPT.imginfo      = 'bdd/image_info/veh_id_stimuli_info.xlsx';
% OPT.outdir       = 'bdd/results/identification/ORIB-data-vb-results-lab-veh-300-rerun/';

%Explanation
% OPT.fixationfile = 'bdd/fixation/veh_exp_screen_10_07_2023_23_49.xlsx';
% OPT.imgdir       = 'bdd/images/orib_veh_id_task_resized/';
% OPT.imgsize      = [1024; 576];
% OPT.imginfo      = 'bdd/image_info/veh_id_stimuli_info.xlsx';
% OPT.outdir       = 'bdd/results/explanation/vehicle/';

% Explanation whole screen
OPT.fixationfile = 'bdd/fixation/veh_exp_screen_10_07_2023_23_49.xlsx';
OPT.imgdir       = 'bdd/images/vehicle_with_screen/';
OPT.imgsize      = [1024; 768];
OPT.imginfo      = 'bdd/image_info/veh_exp_screen_stimuli_info.xlsx';
OPT.outdir       = 'bdd/results/explanation/231015_vehicle_whole_screen_1-group/';
OPT.HEM_K         = 1;    % fit one group to identify outliers


%Categorization and Explanation
% OPT.fixationfile = 'Input-Data/all-cat-exp-data';
% OPT.imgdir       = 'HUAWEI-categorization-data/';
% OPT.imgsize      = [1024; 576];
% OPT.imginfo      = 'HUAWEI-categorization-data/stimuli-info.xlsx';
% OPT.outdir       = 'Categorization-Explanation-results/';
```