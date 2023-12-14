%% This script should be run after obtaining the cogroup hmms for specific
%   edge case stimuli. Set the file name and index of the edge case among 
%   all stimuli. It will replace the corresponding parameters in the
%   cogroup hmms of all cases.

clear

root_dir = ['H:\OneDrive - The University Of Hong Kong\mscoco' '/exp_results_1group_1203/'];
sub_dir = 'edgecase_hem_crghmm_outdir';

%  Original cogrphmms
all_cogroup_filename = 'cogroup_hmms_1203.mat';
all_cogroup_output_filename = sprintf('cogroup_hmms_1203_fixed_%s.mat',datestr(now,'mm_dd_yyyy_HH_MM'));
all_cogroup = load([root_dir all_cogroup_filename]);

edge_filenames = dir(fullfile([root_dir sub_dir],'*.mat')); % load saved coclustering hmms with only one stimulus
for k = 1:length(edge_filenames)
    baseFileName = edge_filenames(k).name;
    fullFileName = fullfile(root_dir, sub_dir, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    
    e_cogroup = load(fullFileName); 

    e = regexp(baseFileName, '\d+','match');
    e = str2num(e{1}); % get sorted index of stimulus in cogrouphmms
    
    % replace corresponding varibles; not sure if all needed

    if contains(fullFileName, 'edgecase_hem')
        % non-vb coclustering 

        all_cogroup.cogroup_hmms{e, 1}.hmms{1, 1} = e_cogroup.cogroup_hmms{1, 1}.hmms{1, 1};
        
        all_cogroup.cogroup_hmms{e, 1}.hemopt.N(e)   = e_cogroup.cogroup_hmms{1, 1}.hemopt.N;
        
    elseif contains(fullFileName, 'edgecase_vbhem')
        all_cogroup.vbco.L_elbo(e,:) = e_cogroup.vbco.L_elbo;
        
        all_cogroup.vbco.cogroup_hmms{1, e} = e_cogroup.vbco.cogroup_hmms{1,1};
        
        all_cogroup.vbco.learn_hyps.vbopt.S(e) = e_cogroup.vbco.learn_hyps.vbopt.S(1);
        
        all_cogroup.vbco.vbhemopt.S(e) = e_cogroup.vbco.vbhemopt.S(1);
    end
end

save([root_dir all_cogroup_output_filename],'-struct','all_cogroup');
