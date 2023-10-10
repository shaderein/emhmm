function [data, sid_names, sid_trials] = read_xls_fixations(xlsname, opt, cropfix, img)
% read_xls_fixations - read an EXCEL file with fixation data
%
%  [data, subject_names, trial_names] = read_xls_fixations(xlsname, opt, cropfix, img)
%
% Expected header cells in the spreadsheet:
%   SubjectID = subject ID
%   TrialID   = trial ID for subject
%   FixX      = fixation X-location
%   FixY      = fixation Y-location
%   FixD      = fixation duration in milliseconds (optional)
%
% SubjectID and TrialID can be either strings or numbers.
% FixX and FixY must be numbers.
% FixD is a number (milliseconds).
%
% Data will be separated by SubjectID and TrialID. 
% For each trial, fixations are assumed to be in sequential order.
%
% INPUT
%   xlsname - filename for the Excel spreedsheet (xls)
%   opt     - options
%                'F' - fixation location only, ignore duration.
%   cropfix - [] - keep all fixations (default)
%             {'r', [x1 y1 x2 y2]} - keep fixations inside of rectangle 
%                                    defined by points (x1,y1) and (x2,y2).
%             {'e', [x y w h]} - keep  fixations inside of an ellipse at (x,y) with width w and height h
%             - can input several regions for keeping, e.g., 
%               {'r', [5,5,20,20], 'r', [55,55,70,70]}
%   img     - image for visualizing the fixation crops [default='']
%
% OUTPUT
%   data - data cell array
%     data{i}         = i-th subject
%     data{i}{j}      = ... j-th trial
%     data{i}{j}(t,:) = ... [x y] location of t-th fixation 
%                        or [x y d] of t-th fixation (location & duration)
%
%   The subject/trials will be assigned numbers according to their order in the spreadsheet.
%   the following two outputs contain the original subject names and trial IDs:
%     subject_names{i} - the subject ID in the spreadsheet for the i-th subject
%     trial_names{i}{j} - the j-th trial ID for i-th subject
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-01-13
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong

% VERSIONS:
%  2017-01-18 - added duration field FixD
%  2017-01-20 - issue an error if FixX, FixY, FixD values are strings.
%  2017-05-25 - fixed bug that converted SID and TID strings into numbers
%  2020-11-29 - v0.78 - added option 'F' - ignore duration
%  2020-11-29 - v0.78 - cropfix - crop fixations, and image visualization
%  2021-07-21 - v0.78 - ignore empty rows

if nargin<2
  opt = '';
end
if nargin<3
  cropfix = {};
end
if nargin<4
  img = '';
end


fprintf('Reading %s\n', xlsname);

% read the XLS file
[ndata, tdata, rdata] = xlsread(xlsname);


% get the headers
headers = {rdata{1,:}};

% find the header indices
SID = find(strcmp('SubjectID', headers));
TID = find(strcmp('TrialID', headers));
FX  = find(strcmp('FixX', headers));
FY  = find(strcmp('FixY', headers));
% v0.78 - no duration option
if any(opt=='F')
  FD = [];
  fprintf('ignoring duration\n');
else
  FD  = find(strcmp('FixD', headers));
end
if length(SID) ~= 1
  error('error with SubjectID');
end
fprintf('- found SubjectID in column %d\n', SID);
if length(TID) ~= 1
  error('error with TrialID');
end
fprintf('- found TrialID in column %d\n', TID);
if length(FX) ~= 1
  error('error with FixX');
end
fprintf('- found FixX in column %d\n', FX);
if length(FY) ~= 1
  error('error with FixY');
end
fprintf('- found FixY in column %d\n', FY);
if length(FD) == 1
  fprintf('- found FixD in column %d\n', FD);
elseif length(FD) > 1
  error('error with FixD -- to many columns');
end

% initialize names and trial names
sid_names = {};
sid_trials = {};
data = {};

% setup the fixation cropping  
ncropfix = 0;
validfix = zeros(0,2);
invalidfix = zeros(0,2);

nempty = 0;

% read data
for i=2:size(rdata,1)
  mysid = rdata{i,SID};
  mytid = rdata{i,TID};
  
  % v0.78: check if empty row
  if isnan(mysid) && isnan(mytid)
    nempty = nempty+1;
    continue
  end  
  
  if isstr(rdata{i,FX})
    error('Value for FixX is text, not a number.');
  end  
  if isstr(rdata{i,FY})
    error('Value for FixY is text, not a number.');
  end
  
  myfxy  = [rdata{i,FX}, rdata{i,FY}];
  
  % v0.78 - crop fixations outside of box
  if ~isempty(cropfix)
    doskip = 1;
    
    % check each crop
    for ff=1:2:length(cropfix)
      switch(cropfix{ff})
        case 'r'
          % rectangle box - check if inside
          bbx = cropfix{ff+1};
          if ( (bbx(1) <= myfxy(1)) && (myfxy(1) <= bbx(3)) && ...
               (bbx(2) <= myfxy(2)) && (myfxy(2) <= bbx(4)) )
            doskip = 0;
          end
    
        case 'e'
          % ellipse - check if inside
          els = cropfix{ff+1};
          if (   (myfxy(1)-els(1))^2/(0.25*els(3)^2) ...
               + (myfxy(2)-els(2))^2/(0.25*els(4)^2) <= 1)            
            doskip = 0;
          end
          
        otherwise
          error('bad crop spec');
      end
    end
    
    % skip this fixation
    if doskip
      ncropfix = ncropfix+1;
      invalidfix(end+1,:) = myfxy;
      continue
    else
      validfix(end+1,:) = myfxy;
    end    
  end
  
  if length(FD) == 1
    % include duration if available    
    if isstr(rdata{i,FD})
      error('Value for FixD is text, not a number.');
    end
    myfxy = [myfxy, rdata{i,FD}];
  end
  
  % convert sid/tid number to string
  if ~ischar(mysid)
    mysid = sprintf('%g', mysid);
  end
  if ~ischar(mytid)
    mytid = sprintf('%g', mytid);
  end
  
  % find subject
  s = find(strcmp(mysid, sid_names));
  if isempty(s)
    % new subject
    sid_names{end+1,1} = mysid;
    sid_trials{end+1,1} = {};
    s = length(sid_names);
    data{s,1} = {};
  end
  
  % find trial
  t = find(strcmp(mytid, sid_trials{s}));
  if isempty(t)
    sid_trials{s,1}{end+1,1} = mytid;
    t = length(sid_trials{s});
    data{s,1}{t,1} = [];
  end
  
  % put fixation
  data{s,1}{t,1}(end+1,:) = myfxy;
end

fprintf('- found %d subjects:\n', length(sid_names));
fprintf('%s ', sid_names{:})
fprintf('\n');
for i=1:length(data)
  fprintf('    * subject %d had %d trials\n', i, length(data{i}));
end

fprintf('- found %d empty rows in the Excel file\n', nempty);

if ~isempty(cropfix)
  fprintf('- ignored %d fixations outside cropfix\n', ncropfix);
 
  figure
  if ~isempty(img)
    im = imread(img);
    imagesc(im);
  end  
  hold on
  for ff=1:2:length(cropfix)
    switch(cropfix{ff})
      case 'r'
        % rectangle
        bbx = cropfix{ff+1};
        plot([bbx(1) bbx(3), bbx(3), bbx(1), bbx(1)], ...
             [bbx(2) bbx(2), bbx(4), bbx(4), bbx(2)], 'g-', 'linewidth', 2);
      case 'e'
        % ellipse
        els = cropfix{ff+1};
        t = linspace(0,2*pi,100);
        xx = 0.5*els(3)*cos(t) + els(1);
        yy = 0.5*els(4)*sin(t) + els(2);
        plot(xx, yy, 'g-', 'linewidth', 2);
        
      otherwise
        error('bad crop spec')
    end
  end
  
  scatter(validfix(:,1), validfix(:,2), 2, 'b.');
  scatter(invalidfix(:,1), invalidfix(:,2), 2, 'r.');  
  hold off
  title(sprintf('kept/cropped fixations (%d/%d)', ...
    size(validfix,1), size(invalidfix,1)));

elseif ~isempty(img)
  % visualize fixations
  figure
  im = imread(img);
  imagesc(im);
  hold on 
  for d=1:length(data)
    tmp = cat(1, data{d}{:});
    scatter(tmp(:,1), tmp(:,2), 2, '.');
  end
  hold off
  title(sprintf('fixations from %s', xlsname));
end




