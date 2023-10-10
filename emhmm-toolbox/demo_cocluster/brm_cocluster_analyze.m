% Data analysis for BRM journal paper experiments
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2021-07-29
% Antoni B. Chan, Janet H. Hsiao, Lan Hui
% City University of Hong Kong, University of Hong Kong

% 2021-07-29: v0.80 - Hui - initial version

clear
close all

if 0
  %% input info (for experiment 10)
  fixationfile = 'brm-data/fixations-new.xlsx';
  imgdir  = 'brm-data/images/';
  imgsize = [640; 280];
  imginfo = 'brm-data/stimuli-info.xlsx';
  outdir  = 'brm-data-vb-results/';
  indhmms = 'individual_hmms.mat';
  grphmms = 'vbcogroup_hmms.mat';

else
  %% input info for exact BRM experiment (experiment 12)
  fixationfile = 'brm-data/fixations-new.xlsx';
  imgdir  = 'brm-data/images/';
  imgsize = [640; 280];
  imginfo = 'brm-data/stimuli-info.xlsx';
  outdir  = 'brm-data-vb-results-paper/';
  indhmms = 'individual_hmms_paper.mat';
  grphmms = 'vbcogroup_hmms_paper.mat';
end


%% BRM original data
% includes all data...except one missing TOL test (computer crashed)
datafile = 'brm-data/metadata_original/AllData-original.xlsx';

% group names
GROUPNAMES = {'Explorative', 'Focused'};

%% load mat file
outfile = [outdir indhmms];
fprintf('loading individual MAT files: %s\n', outfile);
load(outfile);
 
%% load mat file
outfile = [outdir grphmms];
fprintf('load group MAT file: %s\n', outfile);
load(outfile);
HEM_K = vbco.K;


%% column names in Excel file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55%
datacols_name = {'Participant', ...
  'SS_Motivation', 'SR_RT', 'SR_Dprime_Old', 'SR_Dprime_newb', 'SR_Dprime_difference', 'SR_Dprime_overall', ...
  'Nback_V_Dprime', 'Nback_V_ReactionTime', 'Nback_S_Dprime', 'Nback_S_ReactionTime', ...
  'Flanker_Incongruent_ACC','Flanker_Congruent_ACC','Flanker_Neutral_ACC','Flanker_Overall_ACC',...
  'Flanker_Incongruent_RT','Flanker_Congruent_RT','Flanker_Neutral_RT','Flanker_Overall_RT', ...
  'TOL_OutExecutionTime','TOL_OutInitiationTime','TOL_OutMoveCount','TOL_OutMoveScore', ...
  'TOL_OutTimeViolation','TOL_OutTotalCorrectScore','TOL_OutTotalExecutionTime','TOL_OutTotalInitiationTime', ...
  'TOL_OutTotalMoveScore','TOL_OutTotalTimeAll', ...
  'FixAVG', 'FixDurationAVG', 'SacDurationAVG', 'SacDurationTOTAL', 'SacDistanceAVG' ...
  };

datacols_name_ex = {'Participant', ...
  'Pref Rating', 'SR RT', 'SR Dprime Old', 'SR Dprime newb', 'SR Dprime difference', 'SR d'' Overall', ...
  'Nback V Dprime', 'Nback V ReactionTime', 'Nback S Dprime', 'Nback S ReactionTime', ...
  'Flanker Incongruent ACC','Flanker Congruent Acc','Flanker Neutral ACC','Flanker Overall ACC',...
  'Flanker Incongruent RT','Flanker Congruent RT','Flanker Neutral RT','Flanker Overall RT', ...
  'TOL OutExecutionTime','TOL OutInitiationTime','TOL OutMoveCount','TOL OutMoveScore', ...
  'TOL OutTimeViolation','TOL OutTotalCorrectScore','TOL OutTotalExecutionTime','TOL Total Planning Time', ...
  'TOL OutTotalMoveScore','TOL OutTotalTimeAll', ...
  'Avg Num Fix', 'Avg Fix Duration', 'SacDurationAVG', 'Total Sac Duration', 'Avg Sac Length' ...
  };

analysisdir = 'analysis_original/'

%% load data file
[Dnum, Dtxt, Draw] = xlsread(datafile);

figoptsx = {'Position', [100 100 350 300]};
figoptsxx = {'Position', [100 100 200 175]};
figoptsx2 = {'Position', [100 100 270 170]};
figoptsx3 = {'Position', [100 100 500 170]};

% find header column indices
headers = Draw(1,:);
for i=1:length(datacols_name)
  tmp = find(strcmp(datacols_name{i}, headers));
  if isempty(tmp)
    error(sprintf('could not find %s\n', datacols_name{i}));
  end
  datacols_ind(i) = tmp;
end
  
tmpdata = Draw(2:end,datacols_ind);

% extract data
Nsubj = sum(cellfun(@isnumeric, tmpdata(:,1)));
for j=1:length(datacols_name)
  myname = datacols_name{j};
  Mdata.(myname) = zeros(Nsubj,1);
  
  for i=1:Nsubj
    if isnumeric(tmpdata{i,j})
      Mdata.(myname)(i) = tmpdata{i,j};
    else
      Mdata.(myname)(i) = nan;
    end
  end
end


%% make sure that the particpant order matches 
SubjNamesI = cellfun(@str2double, SubjNames);
if any(SubjNamesI ~= Mdata.Participant)
  error('order of subjects do not match');
end

maxlen = max(cellfun(@length, datacols_name));

%% setup output
[success, msg] = mkdir(outdir, analysisdir);

%% do t-test analysis of groups
groupinds = vbco.groups;

outtxt = '=== t-tests ===========================\n';

outtxt = [outtxt sprintf('group sizes:\n%s=%d\n%s=%d\n', GROUPNAMES{1}, length(groupinds{1}), GROUPNAMES{2}, length(groupinds{2}))];

for j=2:length(datacols_name)
    
  myname = datacols_name{j};
  mynamex = datacols_name_ex{j};
  
  for k=1:2
    grp_samples{k} = Mdata.(myname)(groupinds{k});
    grp_mean(k)    = mean(grp_samples{k}, 'omitnan');
    grp_std(k)     = std(grp_samples{k}, 'omitnan');
    grp_len(k)     = sum(isfinite(grp_samples{k}));
    
    grp_stderr(k)  = grp_std(k) / sqrt(grp_len(k));
  end
  
  % run t-test
  [h,p, ci, stats] = ttest2( grp_samples{1},grp_samples{2} );  
  deffect = computeCohen_d(grp_samples{1}, grp_samples{2}, 'independent');
  
  % output text
  if (p<0.05)
    tmp = '** ';
  elseif (p<0.1)
    tmp = '+  ';
  else
    tmp = '   ';
  end
  outtxt = [outtxt tmp];  
  outtxt = [outtxt sprintf('%26s: ', myname)];  
  outtxt = [outtxt sprintf('p=%0.5f; t(%d)=%0.5f; d=%0.5f --- ', p, stats.df, stats.tstat, deffect)];
  for k=1:2
    outtxt = [outtxt sprintf('grp%d(%s): mn=%g; std=%g; ', k, GROUPNAMES{k}, grp_mean(k), grp_std(k))];
  end  
  outtxt = [outtxt '\n'];
  
  % make a plot
  if (p<0.1)
    figure(figoptsxx{:});
    %bar([1], grp_mean(1), 'r');
    plot([1], grp_mean(1), 'ro');
    hold on
    errorbar([1], grp_mean(1), grp_stderr(1), 'r', 'linewidth', 2);
    %bar([2], grp_mean(2), 'b');    
    plot([2], grp_mean(2), 'bo');    
    errorbar([2], grp_mean(2), grp_stderr(2), 'b', 'linewidth', 2);
    hold off
    set(gca, 'XTick', [1 2]);
    set(gca, 'XTickLabel', GROUPNAMES);
    grid on
    ylabel(mynamex);
    set(gca, 'FontSize', 12);
    myax = axis();
    myax(1) = [0.7];
    myax(2) = [2.3];
    axis(myax);
    title(sprintf('p=%0.5f; t(%d)=%0.3f', p, stats.df, stats.tstat));    
    
    % save figure
    outfile = [outdir analysisdir sprintf('ttest_%s', myname)];
    fprintf('saving %s\n', outfile);
    savefigs(outfile, 'png');
  end
end


%% do regression analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% exchange the hmms, just hmms are ok.

% compute log-likelihoods for each subject to each cluster
for i=1:Nsubjects
  for j=1:Nstimuli
    for k=1:HEM_K
      if ~isempty(alldataC{j,i})
        LL(i,j,k) = vbhmm_ll(vbco.cogroup_hmms{j}.hmms{k}, alldataC{j,i});
      else
        LL(i,j,k) = 0;
      end
    end
  end
end
 
% sum over stimuli
LL1 = sum(LL(:,:,1),2);
LL2 = sum(LL(:,:,2),2);

% since some stimuli are missing, then normalize (group 1 is positive,
% group 2 is negative)
HA = (LL1 - LL2) ./ (abs(LL1) + abs(LL2));

outtxt = [outtxt '=== regression ===========================\n'];

%% regression analysis, and make plots
for j=2:length(datacols_name)
    
  myname = datacols_name{j};
  mynamex = datacols_name_ex{j};
  
  [b, bint, r, rint, stats] = regress(Mdata.(myname)(:), [ones(length(HA), 1), HA(:)]);
  R2 = stats(1);
  F  = stats(2);
  p  = stats(3);
  df = sum(isfinite(Mdata.(myname)(:)))-2;
  R = sqrt(R2)*sign(b(2));
  
  [RR,P,RL,RU]  = corrcoef(Mdata.(myname)(:),HA(:));  
  
  if (p<0.05)
    tmp = '** ';
  elseif (p<0.1)
    tmp = '+  ';
  else
    tmp = '   ';
  end
  outtxt = [outtxt tmp];
  
  outtxt = [outtxt sprintf('%26s: ', myname)];
  %outtxt = [outtxt sprintf('p=%0.5f; R=%0.3f; F(%d)=%0.3f', p, R, df, F)];
  outtxt = [outtxt sprintf('p=%0.5f; r(%d)=%0.3f; F=%0.3f', p, df, R, F)];
  outtxt = [outtxt '\n'];
  
  if (p<0.1)
    figure(figoptsx2{:})
    plot(HA(:), Mdata.(myname)(:), 'bo');
    hold on
    xmm = [min(HA); max(HA)];
    plot(xmm, [ones(2,1), xmm]*b, 'k-', 'LineWidth', 2);
    %title(sprintf('p=%0.5f, R=%0.3f; F(%d)=%0.3f', p, R, df, F));    
    title(sprintf('p=%0.5f, r(%d)=%0.3f', p, df, R));
    
    xlabel(sprintf('%s <-- EF Scale --> %s', GROUPNAMES{2}, GROUPNAMES{1}), 'FontWeight', 'Bold');
    %xlabel('Group2 <--         --> Group 1', 'FontWeight', 'Bold');
    
    ylabel(mynamex, 'FontWeight', 'Bold');
    grid on

    outfile = [outdir analysisdir sprintf('regr_%s', myname)];
    fprintf('saving %s\n', outfile);
    savefigs(outfile, 'png');
  end
  
end


%% KL histograms (stimuli type comparison)
grps = vbco.groups;

% sum over subject
KL1 = mean( LL(grps{1},:,1) - LL(grps{1},:,2), 1);
KL2 = mean( LL(grps{2},:,2) - LL(grps{2},:,1), 1);
  
% sort by difference
KLsum = 0.5*(KL1 + KL2);

[KLsort, KLinds] = sort(KLsum(:), 'descend');

meanVehicles = mean(KLsum(1:60));
meanAnimals = mean(KLsum(61:end));
stdVehicles = std(KLsum(1:60));
stdAnimals  = std(KLsum(61:end));

% histogram of all
figure(figoptsx{:})
[HN, HX] = hist(KLsum, 30);
bar(HX, HN)
xlabel('SKL');
ylabel('count');
grid on

outfile = [outdir analysisdir sprintf('stimuli_type_hist')];
fprintf('saving %s\n', outfile);
savefigs(outfile, 'png');

% histogram of each
figure(figoptsx3{:})
[HNv] = hist(KLsum(1:60), HX);
[HNa] = hist(KLsum(61:end), HX);
subplot(2,1,1)
bar(HX, HNv, 'b', 'FaceAlpha', 0.5);
hold on
plot([meanVehicles, meanVehicles], [0 max(HN)+1], 'b--', 'LineWidth', 2);
hold off
grid on
xlabel('SKL')
ylabel('count');
title(sprintf('vehicles (mean=%0.3f, std=%0.3f)', meanVehicles, stdVehicles));
subplot(2,1,2)
bar(HX, HNa, 'r', 'FaceAlpha', 0.5)
hold on
plot([meanAnimals meanAnimals], [0 max(HN)+1], 'r--', 'LineWidth', 2);
hold off
grid on
xlabel('SKL')
ylabel('count');
title(sprintf('animals (mean=%0.3f, std=%0.3f)', meanAnimals, stdAnimals));

outfile = [outdir analysisdir sprintf('stimuli_type_hist2')];
fprintf('saving %s\n', outfile);
savefigs(outfile, 'png');

% scatter plot
figure(figoptsx{:})
plot(1:60, KLsum(1:60), 'bx');
hold on
plot(61:120, KLsum(61:end), 'ro');
hold off
xlabel('image index');
ylabel('SKL');
legend({'vehicles', 'animals'})
grid on

outfile = [outdir analysisdir sprintf('stimuli_type_scatter')];
fprintf('saving %s\n', outfile);
savefigs(outfile, 'png');

[h,p,ci,stats] = ttest2(KLsum(1:60), KLsum(61:end));
deffect = computeCohen_d(KLsum(1:60), KLsum(61:end));

outtxt = [outtxt '=== stimuli test ===================\n'];
outtxt = [outtxt 'vehicles vs animals : SKL between groups \n'];
outtxt = [outtxt sprintf('  p=%0.5f; t(%d)=%0.5f cohen_d=%0.4g --- \n', p, stats.df, stats.tstat, deffect)];
outtxt = [outtxt sprintf('  vehicles: mn=%g; std=%g; ', meanVehicles, stdVehicles)];
outtxt = [outtxt sprintf('  animals:  mn=%g; std=%g; ', meanAnimals, stdAnimals)];
outtxt = [outtxt '\n'];


%% KL comparisons (group comparisons for each stimuli, separately)
stim_h1 = [];
stim_p1 = [];
stim_h2 = [];
stim_p2 = [];
stim_LL1 = [];
stim_LL2 = [];

for j=1:Nstimuli
  % get LL for group 1 under models 1 and 2
  stim_LL1(:,j,1) = LL(groupinds{1}, j, 1);
  stim_LL1(:,j,2) = LL(groupinds{1}, j, 2);
  [stim_h1(j), stim_p1(j)] = ttest2(stim_LL1(:,j,1), stim_LL1(:,j,2), 'tail', 'right');
  
  % get LL for group 2 under models 2 and 1
  stim_LL2(:,j,1) = LL(groupinds{2}, j, 2);
  stim_LL2(:,j,2) = LL(groupinds{2}, j, 1);
  [stim_h2(j), stim_p2(j)] = ttest2(stim_LL2(:,j,1), stim_LL2(:,j,2), 'tail', 'right');  

  % G1 data under G1 model, G2 data under G2 model
  stim_LL1c = [stim_LL1(:,j,1); stim_LL2(:,j,1)];
  % G1 data under G2 model, G2 data under G1 model
  stim_LL2c = [stim_LL1(:,j,2); stim_LL2(:,j,2)];
  [stim_hc(j), stim_pc(j)] = ttest2(stim_LL1c, stim_LL2c, 'tail', 'right');
end

% both
stim_h = stim_h1.*stim_h2;

outtxt = [outtxt '=== group test (each stimuli) ===================\n'];
outtxt = [outtxt 'comparison between group models:\n'];
outtxt = [outtxt sprintf('  number of sig diff models (grp1 data KL): %d/%d \n', sum(stim_h1), length(stim_h1))];
outtxt = [outtxt sprintf('  number of sig diff models (grp2 data KL): %d/%d \n', sum(stim_h2), length(stim_h2))];
outtxt = [outtxt sprintf('  number of sig diff models (both):         %d/%d \n', sum(stim_h), length(stim_h))];
outtxt = [outtxt sprintf('  number of sig diff models (SKL):          %d/%d \n', sum(stim_hc), length(stim_hc))];
  
figure(100)
clf
plot(stim_p1(KLinds), 'rx-');
hold on
plot(stim_p2(KLinds), 'bo-');
plot(0.5*stim_h(KLinds), 'k-');
hold off
xlabel('p');
ylabel('stimuli, sorted by descending SKL')
legend('group1 data', 'group2 data', 'h combined');


%% perform 2-way unbalanced ANOVA
YY = [];
g1 = {};
g2 = {};
for k=1:2
  for j=1:Nstimuli
    for i=1:size(stim_LL1,1)
      YY(end+1) = stim_LL1(i,j,k);
      g1{end+1} = sprintf('n%d', j);
      g2{end+1} = sprintf('g%d', k);
    end
  end
end
[anova1_p, anova1_tbl] = anovan(YY, {g1, g2}, 'model', 2, 'varnames', {'S', 'G'});

YY = [];
g1 = {};
g2 = {};
for k=1:2
  for j=1:Nstimuli
    for i=1:size(stim_LL2,1)
      YY(end+1) = stim_LL2(i,j,k);
      g1{end+1} = sprintf('n%d', j);
      g2{end+1} = sprintf('g%d', k);
    end
  end
end
[anova2_p, anova2_tbl] = anovan(YY, {g1, g2}, 'model', 2, 'varnames', {'S', 'G'});

outtxt = [outtxt '=== group test (all stimuli) ===================\n'];
outtxt = [outtxt 'comparison between group models:\n'];
outtxt = [outtxt 'group 1 anova: \n'];
outtxt = [outtxt print_cell_table(anova1_tbl, 10), '\n\n'];
outtxt = [outtxt 'group 2 anova: \n'];
outtxt = [outtxt print_cell_table(anova2_tbl, 10), '\n\n'];

%% save output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(outtxt);
outfile = [outdir analysisdir 'analysis.txt'];
fprintf('saving %s\n', outfile);
fp = fopen(outfile, 'w');
fprintf(fp, outtxt);
fclose(fp);

%% save data as xlsx %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% put all data into a single structure
Mdata2 = Mdata;
Mdata2.LL1 = LL1;
Mdata2.LL2 = LL2;
Mdata2.HA  = HA;
Mdata2.Group = vbco.label(:);

myT = struct2table(Mdata2);

outfile = [outdir analysisdir 'data_LL.csv'];
fprintf('saving %s\n', outfile);
writetable(myT, outfile)

%% save stimuli KL data
Sdata.index = (1:length(KLsum))';
Sdata.imgname = StimuliNamesC;
Sdata.SKL = KLsum(:);

myT = struct2table(Sdata);
outfile = [outdir analysisdir 'data_KL.csv'];
fprintf('saving %s\n', outfile);
writetable(myT, outfile)

%% sort KL data
Sdata2.index   = Sdata.index(KLinds);
Sdata2.imgname = Sdata.imgname(KLinds);
Sdata2.SKL     = Sdata.SKL(KLinds);

myT = struct2table(Sdata2);
outfile = [outdir analysisdir 'data_KL_sorted.csv'];
fprintf('saving %s\n', outfile);
writetable(myT, outfile)

%% sort p-values
[~, kk] = sort(stim_p1+stim_p2 + 10*(1-stim_h), 'ascend');
Hdata.index = kk(:);
Hdata.imgname = StimuliNamesC(kk);
Hdata.h  = stim_h(kk)';
Hdata.h1 = stim_h1(kk)';
Hdata.p1 = stim_p1(kk)';
Hdata.h2 = stim_h2(kk)';
Hdata.p2 = stim_p2(kk)';

myT = struct2table(Hdata);
outfile = [outdir analysisdir 'data_anova_sorted.csv'];
fprintf('saving %s\n', outfile);
writetable(myT, outfile);


%table_to_Mercury=[];
%for i = 1:Nsubjects
%  for j=1:Nstimuli
%    tmp = (LL(i,j,1) - LL(i,j,2)) ./ (abs(LL(i,j,1)) + abs(LL(i,j,2)));
%    table_to_Mercury=[table_to_Mercury;[LL(i,j,1),LL(i,j,2),tmp]];
%  end
%end
