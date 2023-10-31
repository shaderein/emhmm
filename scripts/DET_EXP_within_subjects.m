% Do it for each class separately: veh vs hum

root = 'H:/OneDrive - The University Of Hong Kong/';

% Pearson Correlation

fprintf("==== Pearson Correlation ====");

DET = load([root 'bdd/results/identification/ORIB-data-vb-results-lab-veh-300-rerun/veh_useful_variables_analysis_300_run.mat']);
EXP = load([root 'bdd/results/explanation/231023_vehicle_posneg_fixed_vhem_alpha/useful_variables_analysis_spilt_half.mat']);

[corr p] = corrcoef(DET.AB, EXP.AB)

[DET_sorted_v, DET_sorted_i] = sort(DET.AB);

mdl = fitlm(DET_sorted_v, EXP.AB(DET_sorted_i));

plot(mdl)
xlabel("DET-AB");
ylabel("EXP-AB");


% Chi-square Test

DET = load([root 'bdd/results/identification/ORIB-data-vb-results-lab-veh-300-rerun/veh_vbcogroup_hmms_300_run.mat']);
EXP = load([root 'bdd/results/explanation/231023_vehicle_posneg_fixed_vhem_alpha/vbcogroup_hmms.mat']);

DET_group1 = DET.cogroup_hmms{1,1}.groups{1};
DET_group2 = DET.cogroup_hmms{1,1}.groups{2};

fprintf("==== Groups in Detection Task ====\nGroup 1: %d Subjects, Group 2: %d Subjects\n",...
    length(DET_group1), length(DET_group2));

EXP_group1 = EXP.vbco.groups{1};
EXP_group2 = EXP.vbco.groups{2};

fprintf("==== Groups in Explanation Task ====\nGroup 1: %d Subjects, Group 2: %d Subjects\n\n",...
    length(EXP_group1), length(EXP_group2));

table = [];

table(1,1) = length(intersect(DET_group1, EXP_group1));
table(1,2) = length(intersect(DET_group1, EXP_group2));
table(2,1) = length(intersect(DET_group2, EXP_group1));
table(2,2) = length(intersect(DET_group2, EXP_group2));

fprintf("==== DET Group 1 ∩ EXP Group 1 ====\n%d common subjects:\n%s\n",...
    length(intersect(DET_group1, EXP_group1)),...
    string(mat2str(intersect(DET_group1, EXP_group1))));

fprintf("==== DET Group 1 ∩ EXP Group 2 ====\n%d common subjects:\n%s\n",...
    length(intersect(DET_group1, EXP_group2)),...
    string(mat2str(intersect(DET_group1, EXP_group2))));

fprintf("==== DET Group 2 ∩ EXP Group 1 ====\n%d common subjects:\n%s\n",...
    length(intersect(DET_group2, EXP_group1)),...
    string(mat2str(intersect(DET_group2, EXP_group1))));

fprintf("==== DET Group 2 ∩ EXP Group 2 ====\n%d common subjects:\n%s\n",...
    length(intersect(DET_group2, EXP_group2)),...
    string(mat2str(intersect(DET_group2, EXP_group2))));


fprintf("==== chi-square ====");

addpath(genpath('../helpers/'));

[h, chi, p] = chi2ind(table, 0.05)