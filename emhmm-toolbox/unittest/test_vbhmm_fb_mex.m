function test_vbhmm_fb_mex(mode)
% test_vbhmm_fb_mex - unit test for vbhmm_fb_mex
%
%   test_vbhmm_fb_mex(mode)
%
%   mode = 1 - unit test [default]
%        = 2 - functional test (compare vbhmm_learn w/ and w/o mex)
%        = 3 - timing test
%
% this is an internal function called after compiling the MEX function
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-08-08
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

if nargin==0
  mode = 1;
end

switch(mode)
  %% unit test
  case 1  
    load test_vbhmm_fb_mex_hmm.mat
  
    opt.usegroups = 0;
    opt.savexi    = 0;
  
    % MEX run
    fbstats_mex = vbhmm_fb(mydata, hmm.varpar, opt);
  
    % MATLAB run
    opt.useMEX = 0;
    fbstats_matlab = vbhmm_fb(mydata, hmm.varpar, opt);
  
    % check consistency
    err = totalerror(fbstats_mex, fbstats_matlab, 4, 0);
    if (err>1e-10)
      error('unit test of vbhmm_fb_mex failed');
    end
    
    %% functional test
  case 2 
    % set random state to be able to replicate results
    rand('state', 101);
    randn('state', 101);
  
    %% Load data from xls %%%%%%%%%%%%%%%%%%%%%%%%%%
    % see the xls file for the format
    [data, SubjNames, TrialNames] = read_xls_fixations('../demo/demodata.xls');
  
    % the data is read and separated by subject and trial, and stored in a cell array:
    % data{i}         = i-th subject
    % data{i}{j}      = ... j-th trial
    % data{i}{j}(t,:) = ... [x y] location of t-th fixation
    
    % the same data is stored in a mat file.
    % load demodata.mat
    
    % the number of subjects
    N = length(data);
      
    %% VB Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    K = 3; % automatically select from K=2 to 3
    vbopt.alpha0 = 0.1;
    vbopt.mu0    = [256;192];
    vbopt.W0     = 0.005;
    vbopt.beta0  = 1;
    vbopt.v0     = 5;
    vbopt.epsilon0 = 0.1;
    vbopt.showplot = 0;
    faceimg = '../demo/face.jpg';
    
    for n=1:10
      mydata = data{n};
      
      rand('state', 101);
      randn('state', 101);
      vbopt.useMEX = 1;
      tic
      hmm = vbhmm_learn(mydata, K, vbopt);
      mex_time = toc
      
      figure(100)
      subplot(2,1,1)
      vbhmm_plot_compact(hmm, faceimg);
      
      % run again with MATLAB implementation
      rand('state', 101);
      randn('state', 101);
      vbopt.useMEX = 0;
      
      tic
      hmm_matlab = vbhmm_learn(mydata, K, vbopt);
      matlab_time = toc
      
      figure(100)
      subplot(2,1,2)
      vbhmm_plot_compact(hmm_matlab, faceimg);
      drawnow
      
      tmphmm = rmfield(hmm, 'vbopt');
      tmphmm_matlab = rmfield(hmm_matlab, 'vbopt');
      err = totalerror(tmphmm, tmphmm_matlab, 1, 0)
      if (err>1e-10)
        error('unit test of vbhmm_fb_mex failed');
      end
      
    end
    
    %% timing test
  case 3    
    load test_vbhmm_fb_mex_hmm.mat  
    opt.usegroups = 0;
    opt.savexi    = 0;
  
    %% random data
    NTrials = 1000;
    
    all_data = cell(1,NTrials);
    all_fbstats_mex = cell(1,NTrials);
    all_fbstats_matlab = cell(1,NTrials);
    
    for n=1:NTrials
      all_data{n} = data{mod(n-1,10)+1};
      for i=1:length(all_data{n})
        all_data{n}{i} = all_data{n}{i}+randn(size(all_data{n}{i}));
      end
    end
    
    %% timing run
    % mex (array indexing) = 0.4548
    % matlab = 15.4227
    tic
    for trial=1:NTrials
      % MEX run
      opt.useMEX = 1;
      all_fbstats_mex{trial} = vbhmm_fb(all_data{trial}, hmm.varpar, opt);
    end
    mextime = toc
    
    tic
    for trial=1:NTrials
      % MATLAB run
      opt.useMEX = 0;
      all_fbstats_matlab{trial} = vbhmm_fb(all_data{trial}, hmm.varpar, opt);
    end
    matlabtime = toc
    
    for n=1:NTrials
      tmp = totalerror(all_fbstats_mex{n}, all_fbstats_matlab{n}, 4, 0);
      if (tmp>1e-10)
        error('unit test of vbhmm_fb_mex failed');
      end
    end
end

