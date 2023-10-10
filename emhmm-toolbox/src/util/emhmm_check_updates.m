function emhmm_check_updates()
% emhmm_check_updates - check for updates to the toolbox
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-05-25
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong

% 2017-06-06: ABC - check if this version number is after web version

% useful URLs
url_version   = 'http://visal.cs.cityu.edu.hk/static/downloads/emhmm/VERSION.txt';
url_changelog = 'http://visal.cs.cityu.edu.hk/static/downloads/emhmm/CHANGELOG.txt';
url_project   = 'http://visal.cs.cityu.edu.hk/research/emhmm/';

% check current version
v = emhmm_version();
fprintf('*** <strong>EMHMM - Eye-Movement analysis with HMMs</strong> ***\n');
fprintf('  - current version: %s\n', v);

fprintf('  - checking for new version...\n');

% get version on server
newv = deblank(get_content(url_version));
if isempty(newv);
  fprintf('  - could not contact server, or no internet connection.\n');

else
  
  % check with current version
  if strcmp(newv, v)
    % same version
    fprintf('  - version up-to-date.\n');
    
  else
    % get the actual numbers
    newv_num = str2num(newv(2:end));
    v_num    = str2num(v(2:end));
    
    % this version is more recent than the official version
    if (newv_num < v_num)
      fprintf('  - official version is %s\n', newv);
      fprintf('  - you are using a preview version %s\n', v);
      
    % our version is older than the web version
    else
      fprintf('*** <strong>a new version %s is available!</strong> ***\n', newv);
      
      % get change log
      clog = deblank(get_content(url_changelog));
      if isempty(clog)
        fprintf('  - could not contact server, or no internet connection.\n');
      else
        newinfo = parse_log(clog, v);
        
        fprintf('*** Updates since your version ***\n');
        
        for i=1:length(newinfo)
          fprintf(['  ' newinfo{i} '\n']);
        end
      end
      
      % message for downloading
      fprintf('*** Please visit the project website to download the new version ***\n');
      fprintf('  <a href="%s">%s</a>\n', url_project, url_project);
    end
  end
  
end


%% get content from URL %%%
function out = get_content(url)
if exist('webread')
  try 
    out = webread(url);
  catch me
    out = '';
  end
  
elseif exist('urlread')
  try
    out = urlread(url);
  catch me
    out = '';
  end
  
else
  warning('Your version of MATLAB is way too old!');
end

%% parse CHANGELOG %%%
function out = parse_log(log, curver)

% parse each line
slog = strsplit(log, '\n');
N = length(slog);
out = {};

% find changelog line
i = 1;
while (i<=N)
  if strcmp(slog{i}, '--- CHANGE LOG ---')
    i = i+1;
    break;
  end
  i = i+1;
end

% loop until our version appears
while (i<=N)
  myline = slog{i};
  
  % start of a version number
  if (myline(1) == 'v')
    tmp = strsplit(myline, ' ');
    if strcmp(tmp{1}, curver)
      % this is the current version, so break the loop.
      break;
    end 
  end
  
  % add this line to output
  out{end+1} = myline;
  
  i = i+1;
end

