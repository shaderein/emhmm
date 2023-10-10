% setup the path for the toolbox, and check for updates
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-01-13
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong

% 2016-05-25: ABC - added check for updates
% 2017-08-02: ABC - added check for MEX files


% add path using full path names
myname = mfilename('fullpath');

[pathstr,name,ext] = fileparts(myname);

gdir = [pathstr filesep 'src'];
ddir = [pathstr filesep 'demo'];

addpath(genpath(gdir))

% check for updates
emhmm_check_updates()

% check for MEX files
emhmm_mex_check()

% set global variables

% set default of color mode
global EMHMM_COLORPLOTS
EMHMM_COLORPLOTS = 1;