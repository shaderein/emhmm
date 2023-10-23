function remove_outliers(alldataC, cg_hmms, SubjNames, StimuliNamesImagesC, TrialNamesC, OPT)

% To use the figmerge helper function
addpath(genpath('../helpers/'));

% input data
%   - data - data cell array
%       data{j}         = j-th trial (stimulus)
%       data{j}{i}      = ... i-th subject
%       data{j}{i}(t,:) = ... [x y] location of t-th fixation 
%                          or [x y d] of t-th fixation (location & duration)
%   - cg_hmms - co-clustering results
%
%   Mapping between data idx and sub/stimuli/trial ids
%   - SubjNames
%   - StimuliNamesImagesC
%   - TrialNamesC
%
%   - OPT - for files saving

% output data
% cleaned fixation reports
% dimension

% coordinates of the region within which we will find the center of 
%   ROI to remove
% top-left, top-right, bottom-right, bottom-left
target_ROI = [0,500; 230,500; 230,670; 0,670];

class = 'veh';

figopts = {'Position', [20 0 1240 800], 'visible','off'};

found_ROI_dir = '7_found_textbox_ROI/';
[success, msg] = mkdir(OPT.outdir, found_ROI_dir);

alldataC_cleaned = alldataC;

for j=1:length(StimuliNamesImagesC)

    group_hmms = cg_hmms{j};
    data = alldataC(j,:);
    img = StimuliNamesImagesC{j,2};
    img_num = str2num(erase(StimuliNamesImagesC{j,1},".jpg"));

    K = length(group_hmms.hmms);
    S = length(group_hmms.hmms{1}.prior);
    D = length(group_hmms.hmms{1}.pdf{1}.mean);

    % read gphmm
    % identify ROI for text box (by mean's position?)
    targets = [];
    for r=1:S
        ROI = group_hmms.hmms{1}.pdf{r}.mean;

        if target_ROI(1,1) < ROI(1) & ROI(1) < target_ROI(2,1) &... % x
            target_ROI(1,2) < ROI(2) & ROI(2) < target_ROI(3,2)  % y
            targets = [targets, r];
        end
    end
    if length(targets) > 1
        fprintf('Stimulus %d has more than 1 text-box ROI detected!\n',img_num);
    elseif isempty(targets)
        fprintf('Stimulus %d has 0 text-box ROI detected!\n',img_num);
        
        hfig_orig = vhem_plot_fixations(data, group_hmms, img, 'c', 1, figopts);
        h = findobj(gcf,'type','axes');
        orig_title = get(get(h,'title'),'string');
        set(get(h,'title'),"string",sprintf('(Before) Class=%s Image=%d No ROI',class, img_num));
         
        outfile = [OPT.outdir found_ROI_dir sprintf('Stimuli_%04d', img_num)];
        % fprintf('saving co-clustering group HMMs: %s\n', outfile);
        export_fig(outfile, '-png'); %'-transparent'
        continue
    end
    target = targets(1); % text-box ROI found

    % Find cluster that each fixation belogs to
    % Adapted from vhem_plot_fixations.m
    gamma = cell(1,length(data));

    % for each subject
    for i=1:length(data)
      mydata = data{i};  

      % 2018-11-24: v0.74 - handle empty data
      if isempty(mydata)
        mygamma = {};
        ss = {};
      else
        k = group_hmms.label(i);

        % get most-likely ROI sequences
        myhmm = group_hmms.hmms{k};
        ss = vbhmm_map_state(myhmm, mydata);

        % setup as gamma for plot_emissions
        mygamma = cell(1,length(mydata));
        for jj=1:length(mydata)
          mygamma{jj} = zeros(S, size(mydata{jj},1));
          for q=1:size(mydata{jj},1)
            mygamma{jj}(ss{jj}(q), q) = 1;
          end
        end
      end

      gamma{i} = mygamma; % mapping from fixations to cluster

      % find column (fixation) clustered into target ROI. delete
      % don't remove if only two ROIs exists. (may incorrectly remove the
      % explorative ROI that overlaps with the text-box ROI)
    
      if strcmp(class,'hum') && img_num == 1578
         % special case: hum_1578.jpg has two text-box ROI 3&4
        removed = [find(mygamma{1}(3,:)==1) find(mygamma{1}(4,:)==1)];
        alldataC_cleaned{j,i}{1}(removed,:) = [];
      elseif strcmp(class,'hum') && (img_num == 1066 || img_num == 1624)
         % special case: hum_1066.jpg has 1 text-box ROI2 and 1 actual ROI
        removed = find(mygamma{1}(2,:)==1);
        alldataC_cleaned{j,i}{1}(removed,:) = [];
      elseif ~isempty(mydata) && S > 2
        removed = find(mygamma{1}(target,:)==1);
        alldataC_cleaned{j,i}{1}(removed,:) = [];
      end 

      % Edge case: all fixations fall in this ROI and were deleted
      %   (i.e. the subject only focused on the textbox when viewing)
      % Delete the cell with empty array. Otherwise will cause idx out
      % of bound error when plotting cleaned fixations
      if ~isempty(mydata) && length(alldataC_cleaned{j,i}{1})==0
          alldataC_cleaned{j,i} = [];
      end

    end

    % plot clusters with fixations before & after removal

    % before outliers removal
    hfig_orig = vhem_plot_fixations(data, group_hmms, img, 'c', 1, figopts);
    h = findobj(gcf,'type','axes');
    orig_title = get(get(h,'title'),'string');
    set(get(h,'title'),"string",sprintf('(Before) Class=%s Image=%d ROI=%s',class, img_num, sprintf('%d ',targets)));

    % after outliers removal
    hfig = vhem_plot_fixations(alldataC_cleaned(j,:), group_hmms, img, 'c', 1, figopts);
    h = findobj(gcf,'type','axes');
    orig_title = get(get(h,'title'),'string');
    set(get(h,'title'),"string",sprintf('(Removed) Class=%s Image=%d ROI=%s',class, img_num, sprintf('%d ',targets)));

    figmerge([hfig_orig,hfig],[2,0],1,false,false);

    %drawnow
    outfile = [OPT.outdir found_ROI_dir sprintf('Stimuli_%04d_ROI_%d', img_num, targets(1))];
    % fprintf('saving co-clustering group HMMs: %s\n', outfile);
    export_fig(outfile, '-png'); %'-transparent'
end

% convert data back into fixation report format
[count_img, count_subject] = size(alldataC_cleaned);
for j=1:count_img
    for i=1:count_subject
        if isempty(alldataC_cleaned{j,i})
            continue
        end
        [count_fix,~] = size(alldataC_cleaned{j,i}{1});
        alldataC_cleaned{j,i} = cat(2,...
                                    repmat(convertCharsToStrings(SubjNames{i}),[count_fix,1]),... %SubjectID
                                    repmat(TrialNamesC{j,i},[count_fix,1]),...%TrialID
                                    repmat(StimuliNamesImagesC{j,1},[count_fix,1]),...%StimuliID
                                    num2cell(alldataC_cleaned{j,i}{1})); %fixation X Y
    end
end

% ignore empty cells (subject i doesn't have fixation in this trial OR
%   fixation outside the screen
idx = cellfun('isempty',alldataC_cleaned);
flatten = alldataC_cleaned(~idx);
flatten = vertcat(flatten{:});

flatten = array2table(flatten, 'VariableNames',{'SubjectID','TrialID','StimuliID','FixX','FixY'});
% convert fixation coordinates into number from strings
flatten.FixX = double(flatten.FixX);
flatten.FixY = double(flatten.FixY) - 30; % Realigned fixations with image

filename = sprintf('%s_exp_cleaned_outliers_%s.xlsx', class, datestr(now,'mm_dd_yyyy_HH_MM'));
writetable(flatten, [OPT.outdir '/' filename]);