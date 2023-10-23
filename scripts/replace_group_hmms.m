%% This script should be run after obtaining the cogroup hmms for specific
%   edge case stimuli. Set the file name and index of the edge case among 
%   all stimuli. It will replace the corresponding parameters in the
%   cogroup hmms of all cases.

clear

% vehicle
dir = '../bdd/results/explanation/231018_vehicle_whole_screen_vb_fixed_pos/';
edge_filenames = ["vbcogroup_hmms_71_329jpg_fixed_HEM_S_4.mat", ...
                    "vbcogroup_hmms_102_570jpg_fixed_HEM_S_4.mat", ...
                    "vbcogroup_hmms_113_67jpg_fixed_HEM_S_3.mat"];
edge_ind = [71,102,113];

% human

%%
edge_cases = dictionary(edge_ind, edge_filenames);

all_cogroup_filename = 'vbcogroup_hmms_orig.mat';
all_cogroup_output_filename = sprintf('vbcogroup_hmms_pos_fixed_%s.mat',datestr(now,'mm_dd_yyyy_HH_MM'));
all_cogroup = load([dir all_cogroup_filename]);

all_cogroup_orig = all_cogroup; % debug

%%

for e = keys(edge_cases)'
    e_filename = edge_cases(e);
    e_cogroup = load(strcat(dir,e_filename));

    % replace corresponding varibles; not sure if all needed

    all_cogroup.vbco.L_elbo(e,:) = e_cogroup.vbco.L_elbo;

    all_cogroup.vbco.cogroup_hmms{1, e} = e_cogroup.vbco.cogroup_hmms{1,1};

    all_cogroup.vbco.learn_hyps.vbopt.S(e) = e_cogroup.vbco.learn_hyps.vbopt.S(1);

    all_cogroup.vbco.vbhemopt.S(e) = e_cogroup.vbco.vbhemopt.S(1);
end

save([dir all_cogroup_output_filename],'-struct','all_cogroup');
