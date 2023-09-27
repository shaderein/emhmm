function [outdata] = preprocess_fixations(data, StimuliNames, bboxfile, StimuliNamesImages)
% preprocess_fixations - remove fixations outside an image or bounding box.
%
%  This is used for co-clustering, to remove fixations outside of the images.
%
%  [outdata] = preprocess_fixations(data, StimuliNames, bboxfile, StimuliNamesImages)
%
%  INPUTS
%                data - raw fixation data
%        StimuliNames - names of the stimuli
%            bboxfile - name of the excel file with image sizes (see demo_cocluster for example)
%   StkmuliNameImages - the actual images of each stimuli
%
%  OUTPUTS
%             outdata - the pre-processed fixations
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2019-07-13
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

% 2020-03-12: v0.77 - initial version


fprintf('Preprocessing fixations using %s:\n', bboxfile)

%% read xls file
[bbox_nums, bbox_txt, bbox_raw] = xlsread(bboxfile);

% get headers
STID = find(strcmp('StimuliID', bbox_raw(1,:)));
XLO  = find(strcmp('Xlo', bbox_raw(1,:)));
XHI  = find(strcmp('Xhi', bbox_raw(1,:)));
YLO  = find(strcmp('Ylo', bbox_raw(1,:)));
YHI  = find(strcmp('Yhi', bbox_raw(1,:)));

% read contents
box_xlo = cell2mat(bbox_raw(2:end,XLO));
box_xhi = cell2mat(bbox_raw(2:end,XHI));
box_ylo = cell2mat(bbox_raw(2:end,YLO));
box_yhi = cell2mat(bbox_raw(2:end,YHI));
STnames = bbox_raw(2:end,STID);

outdata = cell(size(data));

h = figure;

% for each subject
for i=1:length(data)
  mydata = data{i};
  mynewdata = cell(size(mydata));
  
  % for each trial
  for j=1:length(mydata)
    % get stimuli name
    mystim = StimuliNames{i}{j};
    k = find(strcmp(mystim, STnames));
    if isempty(k)
      error(sprintf('unknown Stimuli %s', mystim));
    end
    
    % get box
    xlo = box_xlo(k);
    xhi = box_xhi(k);
    ylo = box_ylo(k);
    yhi = box_yhi(k);
    
    % check boundaries
    mytrial = mydata{j};
    S =   (xlo <= mytrial(:,1)) & (mytrial(:,1) <= xhi) ...
        & (ylo <= mytrial(:,2)) & (mytrial(:,2) <= yhi);
    
    ii = find(~S);
    if ~isempty(ii)     
      % show the bad fixations
      if ~isempty(StimuliNamesImages)
        figure(h)
        clf
        kk = find(strcmp(mystim, StimuliNamesImages(:,1)));
        showim(StimuliNamesImages{kk,2});
        hold on
        plot([xlo xhi xhi xlo xlo], [ylo ylo yhi yhi ylo], 'g-');
        plot(mytrial(S,1), mytrial(S,2), 'bo');
        plot(mytrial(ii,1), mytrial(ii,2), 'rx');
        hold off
        drawnow
      end
      
      % remove fixation
      fprintf('- Subject %d Trial %d: removing %d fixations (%d remaining)\n', i, j, length(ii), sum(S));
      if sum(S) > 0
        mytrial = mytrial(S,:);
      else
        mytrial = [];
      end
      
      
      %if length(mytrial)==0
      %  error('no fixations for this trial!');
      %end
      
    end
    
    % store the trial
    mynewdata{j} = mytrial;       
  end
  
  % store the data
  outdata{i} = mynewdata;
end