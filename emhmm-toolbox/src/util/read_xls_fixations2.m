function [data, sid_names, sid_trials, sid_stimuli] = read_xls_fixations2(xlsname, opt)
% read_xls_fixations2 - read an EXCEL file with fixation data  (for co-clustering)
%
%  [data, subject_names, trial_names, stimuli_names] = read_xls_fixations(xlsname, opt)
%
% Expected header cells in the spreadsheet:
%   SubjectID = subject ID
%   TrialID   = trial ID for subject
%   FixX      = fixation X-location
%   FixY      = fixation Y-location
%   FixD      = fixation duration in milliseconds (optional)
%   StimuliID = stimuli name (optional)
%
% SubjectID and TrialID can be either strings or numbers.
% FixX and FixY must be numbers.
% FixD is a number (milliseconds).
% StimuiID can be either strings or numbers.
%
% Data will be separated by SubjectID and TrialID. 
% For each trial, fixations are assumed to be in sequential order.
%
% INPUT
%   xlsname - filename for the Excel spreedsheet (xls)
%   opt     - options (not used)
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
%     stimuli_names{i}{j} - the stimuliID for the j-th trial and i-th subject
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2019-07-13
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

% VERSIONS:
%  2017-01-18 - added duration field FixD
%  2017-01-20 - issue an error if FixX, FixY, FixD values are strings.
%  2017-05-25 - fixed bug that converted SID and TID strings into numbers
%  2018-06-14 - v0.74 - add support for StimuliID in excel file.
%  2020-03-12: v0.77 - initial version
%  2021-05-19: v0.80 - skip empty rows (NaN for SubjectID), issue error of NaN for fixation values.

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
FD  = find(strcmp('FixD', headers));
STID = find(strcmp('StimuliID', headers));

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
  error('error with FixD -- too many columns');
end
if length(STID) == 1
  fprintf('- found StimuliID in column %d\n', STID);
end

% initialize names and trial names
sid_names = {};
sid_trials = {};
data = {};
sid_stimuli = {};

nempty = 0;

% read data
for i=2:size(rdata,1)
  mysid = rdata{i,SID};
  mytid = rdata{i,TID};
  
  % v0.XX: check if empty row
  if any(isnan(mysid)) && any(isnan(mytid))
    nempty = nempty+1;
    continue
  end
  
  if isstr(rdata{i,FX})
    error('Value for FixX is text, not a number.');
  end  
  if isstr(rdata{i,FY})
    error('Value for FixY is text, not a number.');
  end
  
  if isnan(rdata{i,FX})
    error('Value for FixX is NaN')
  end
  if isnan(rdata{i,FY})
    error('Value for FixY is NaN')
  end
  
  myfxy  = [rdata{i,FX}, rdata{i,FY}];
  if length(FD) == 1
    % include duration if available    
    if isstr(rdata{i,FD})
      error('Value for FixD is text, not a number.');
    end
    myfxy = [myfxy, rdata{i,FD}];
  end
  
  if length(STID) == 1
    mystid = rdata{i, STID};
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
    if length(STID)==1
      sid_stimuli{end+1,1} = {};
    end
  end
  
  % find trial
  t = find(strcmp(mytid, sid_trials{s}));
  if isempty(t)
    % add new trial
    sid_trials{s,1}{end+1,1} = mytid;
    
    % start fixation sequence
    t = length(sid_trials{s});
    data{s,1}{t,1} = [];
    
    % add stimuli for this trial
    if length(STID)==1
      sid_stimuli{s,1}{end+1,1} = mystid;
    end
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


