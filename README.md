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