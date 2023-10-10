function emhmm_mex_check(force)
% emhmm_mex_check - check MEX files in toolbox
%
%  emhmm_mex_check(force)
%
%    force = [] - skip compiling if MEX files exist [default]
%            'all' - force compile and unit test
%            'funcname' - force compile of a specific function 'funcname'
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-08-02
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong

% v0.80 - ABC - added VBHEM mex file

if (nargin==0)
  force = [];
end

% change to the base directory
OLDDIR = pwd;
cd(emhmm_base_dir());

try
  
  % list of MEX files and locations
  emhmm_mex(1).func = 'vbhmm_fb_mex';
  emhmm_mex(1).file = 'vbhmm_fb_mex.c';
  emhmm_mex(1).dir  = 'src/hmm';
  emhmm_mex(1).unittest = 'test_vbhmm_fb_mex';
  
  emhmm_mex(2).func = 'hem_hmm_bwd_fwd_mex';
  emhmm_mex(2).file = 'hem_hmm_bwd_fwd_mex.c';
  emhmm_mex(2).dir  = 'src/hem/vhem_h3m';
  emhmm_mex(2).unittest = 'test_hem_hmm_bwd_fwd_mex';
    
  emhmm_mex(3).func = 'vbhem_hmm_bwd_fwd_mex';
  emhmm_mex(3).file = 'vbhem_hmm_bwd_fwd_mex.c';
  emhmm_mex(3).dir  = 'src/vbhem';
  emhmm_mex(3).unittest = 'test_vbhem_hmm_bwd_fwd_mex';
  
  fprintf('  - checking MEX files...\n');
  
  allokay = 1;
  
  for i=1:length(emhmm_mex)
    myfunc = emhmm_mex(i).func;
    mexhere = (exist(myfunc)==3);
    
    if mexhere && (isempty(force) || (~strcmp(force, 'all') && ~strcmp(force, myfunc)))
      % found mex ... do nothing
      %fprintf('  * %s found\n', myfunc);
      
    else
      allokay = 0;
      
      if (~mexhere)
        % didn't find mex file
        fprintf('    * %s not found\n', emhmm_mex(i).func);
      end
      if (~isempty(force))
        fprintf('    * force compile of %s\n', emhmm_mex(i).func);
      end
      
      % change to directory of MEX file
      thisdir = pwd;
      cd(emhmm_mex(i).dir);
      
      % try to compile
      try
        fprintf('      + compiling...\n--\n');
        mex(emhmm_mex(i).file);
        fprintf('--\n      + done\n');
        hadproblem = 0;
      catch me
        fprintf('       + problem compiling MEX file:\n--\n');
        me.message
        fprintf('--\n       + falling back to MATLAB implementation. This will be slower.\n');
        hadproblem = 1;
      end
      
      % go back to where we started
      cd(thisdir);
      
      % unit test
      if (~hadproblem)
        cd('unittest');
        feval(emhmm_mex(i).unittest);
        fprintf('      + unit test OK\n');
      end
      
      % go back to where we started
      cd(thisdir);
    end
  end
  
  if(allokay)
    fprintf('  - all OK\n');
  end

  % catch errors
catch me
  
  % go back to original directory
  cd(OLDDIR);
  
  % re-throw the error
  me.rethrow();
end

% go back to original directory
cd(OLDDIR);
