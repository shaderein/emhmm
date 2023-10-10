function [basedir] = emhmm_base_dir()
% emhmm_base_dir - get the base directory of EMHMM
%
%  basedir = emhmm_base_dir()
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-08-17
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

myname = mfilename('fullpath');

% this file should be:
%  <basedir>/src/util/emhmm_basedir.m

[pathstr,name,ext] = fileparts(myname);

[pathstr2,name2,ext] = fileparts(pathstr);
if ~strcmp(name2, 'util')
  error('emhmm_base_dir is misplaced?');
end

[pathstr3,name3,ext] = fileparts(pathstr2);
if ~strcmp(name3, 'src')
  error('emhmm_base_dir is misplaced?');
end
  

basedir = pathstr3;