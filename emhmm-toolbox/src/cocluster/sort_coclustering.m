function [alldataC, TrialNamesC, StimuliNamesC] = sort_coclustering(alldata, TrialNames, StimuliNames)
% sort_coclustering - sort data for co-clustering
%
%  [alldataC, TrialNamesC, StimuliNamesC] = sort_coclustering(alldata, TrialNames, StimuliNames)
%
%  sort data by collecting by stimuli name
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2019-07-13
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

% 2020-03-12: v0.77 - initial version


% handle missing stimuli

% get the list of stimuli from all subjects
tmp = {};
for i=1:length(StimuliNames)
  tmp = [tmp; StimuliNames{i}(:)];
end

StimuliNamesC = unique(tmp);
Nsubjects = length(alldata);
Nstimuli  = length(StimuliNamesC);


fprintf('Sorting for co-clustering:\n');
fprintf('NOTE: stimuli order is changed\n');
fprintf('- Found %d stimuli: ', Nstimuli);
for i=1:Nstimuli
  fprintf('%s ', StimuliNamesC{i});
end
fprintf('\n');

% new cell array: stimuli x subjects
alldataC = cell(Nstimuli, Nsubjects);
TrialNamesC = cell(Nstimuli, Nsubjects);

% go through each subject
for i=1:Nsubjects
  mydata = alldata{i};
  
  % go through each data
  for j=1:length(mydata)
    
    % find the stimuli in the master list
    mystim = StimuliNames{i}{j};
    stimind = find(strcmp(mystim, StimuliNamesC));
    
    if isempty(stimind)
      error(sprintf('subject %d has unknown stimuli "%s"', i, mystim));
    end
    
    if ~isempty(mydata{j})
      % store data and trialname
      alldataC{stimind,i}{end+1} = mydata{j};
      TrialNamesC{stimind,i}{end+1} = TrialNames{i}{j};
    end
  end
end

% check if all stimuli are present
tmpstr = '';
for i=1:Nsubjects
  fprintf('  * subject %d trials: ', i);
  for j=1:Nstimuli
    tmp = length(alldataC{j,i});
    if (tmp==0)
      tmpstr = [tmpstr, sprintf('- subject %d is missing stimuli %s (ind=%d)\n', i, StimuliNamesC{j}, j)];
    end
    fprintf('%3d ', tmp);
  end
  fprintf('\n');
end
if length(tmpstr)>0
  fprintf(tmpstr);
end
