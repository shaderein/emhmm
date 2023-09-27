function test_hem_hmm_bwd_fwd_mex(mode)
% test_hem_hmm_bwd_fwd_mex - unit test for hem_hmm_bwd_fwd_mex
%
%   test_hem_hmm_bwd_fwd_mex(mode)
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
if (nargin==0)
  mode = 1;
end

switch(mode)
  case 1
    %% unit test on MEX function (diagonal covariance)
    load test_hem_hmm_bwd_fwd_mex_h3ms.mat
        
    %% MEX function
    [xL_elbo, xsum_nu_1, xupdate_emit_pr, xupdate_emit_mu, xupdate_emit_Mu, xsum_xi] = ...
      hem_hmm_bwd_fwd_mex(h3m_b.hmm,h3m_r.hmm,T,smooth, 5, 5);
    
    %% MATLAB function
    % storage
    L_elbo           = zeros(Kb,Kr);
    sum_nu_1         = cell(Kb,Kr);
    update_emit_pr   = cell(Kb,Kr);
    update_emit_mu   = cell(Kb,Kr);
    update_emit_Mu   = cell(Kb,Kr);
    sum_xi           = cell(Kb,Kr);
    
    % loop reduced mixture components
    for j = 1 : Kr
      % get this HMM in the reduced mixture
      hmm_r = h3m_r.hmm{j};
      % loop base mixture components
      for i = 1 : Kb
        % get this HMM in the base mixture
        hmm_b = h3m_b.hmm{i};
        
        [L_elbo(i,j) ...
          sum_nu_1{i,j} ...
          update_emit_pr{i,j} ...
          update_emit_mu{i,j} ...
          update_emit_Mu{i,j}  ...
          sum_xi{i,j} ...
          ] = hem_hmm_bwd_fwd(hmm_b,hmm_r,T,smooth, h3m_b_hmm_emit_cache{i});
      end
    end
    
    % consistency check
    err = [];
    err(end+1) = totalerror(xL_elbo, L_elbo, 4, 0);
    err(end+1) = totalerror(xsum_nu_1, sum_nu_1, 4, 0);
    err(end+1) = totalerror(xsum_xi, sum_xi, 4, 0);
    err(end+1) = totalerror(xupdate_emit_pr, update_emit_pr, 4, 0);
    err(end+1) = totalerror(xupdate_emit_mu, update_emit_mu, 4, 0);
    err(end+1) = totalerror(xupdate_emit_Mu, update_emit_Mu, 4, 0);
    if (any(err>1e-10))
      err
      error('unit test for test_hem_hmm_bwd_fwd_mex failed');
    end
  
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% unit test on MEX function (full covariance)
    clear
    load test_hem_hmm_bwd_fwd_mex_h3ms_full.mat
        
    %% MEX function
    [xL_elbo, xsum_nu_1, xupdate_emit_pr, xupdate_emit_mu, xupdate_emit_Mu, xsum_xi] = ...
      hem_hmm_bwd_fwd_mex(h3m_b.hmm,h3m_r.hmm,T,smooth, maxN, maxN2, logdetCovR, invCovR);
    
    %% MATLAB function
    % storage
    L_elbo           = zeros(Kb,Kr);
    sum_nu_1         = cell(Kb,Kr);
    update_emit_pr   = cell(Kb,Kr);
    update_emit_mu   = cell(Kb,Kr);
    update_emit_Mu   = cell(Kb,Kr);
    sum_xi           = cell(Kb,Kr);
    
    % loop reduced mixture components
    for j = 1 : Kr
      % get this HMM in the reduced mixture
      hmm_r = h3m_r.hmm{j};
      % loop base mixture components
      for i = 1 : Kb
        % get this HMM in the base mixture
        hmm_b = h3m_b.hmm{i};
        
        [L_elbo(i,j) ...
          sum_nu_1{i,j} ...
          update_emit_pr{i,j} ...
          update_emit_mu{i,j} ...
          update_emit_Mu{i,j}  ...
          sum_xi{i,j} ...
          ] = hem_hmm_bwd_fwd(hmm_b,hmm_r,T,smooth, h3m_b_hmm_emit_cache{i});
      end
    end
    
    % consistency check
    err = [];
    err(end+1) = totalerror(xL_elbo, L_elbo, 4, 0);
    err(end+1) = totalerror(xsum_nu_1, sum_nu_1, 4, 0);
    err(end+1) = totalerror(xsum_xi, sum_xi, 4, 0);
    err(end+1) = totalerror(xupdate_emit_pr, update_emit_pr, 4, 0);
    err(end+1) = totalerror(xupdate_emit_mu, update_emit_mu, 4, 0);
    err(end+1) = totalerror(xupdate_emit_Mu, update_emit_Mu, 4, 0);
    if (any(err>1e-10))
      err      
      error('unit test for test_hem_hmm_bwd_fwd_mex failed');
    end
    
    
  case 2
    %% functional test (diagonal)
    load test_hem_hmm_bwd_fwd_mex_hmms.mat
  
    hemopt.sortclusters = 'f';
    hemopt.trials = 5;
    hemopt.emit.covar_type = 'diag';
    faceimg = '../demo/face.jpg';
    
    % cluster
    rand('state', 101);
    randn('state', 101);
    
    hemopt
    tic
    [xgroup_hmms2] = vhem_cluster(hmm, 2, 3, hemopt);
    vhem_plot(xgroup_hmms2, faceimg);
    mex_time = toc
  
    rand('state', 101);
    randn('state', 101);
    tic
    hemopt.useMEX = 0
    [group_hmms2] = vhem_cluster(hmm, 2, 3, hemopt);
    vhem_plot(group_hmms2, faceimg);
    matlab_time = toc
  
    tmpgroup        = rmfield(xgroup_hmms2, 'hemopt');
    tmpgroup_matlab = rmfield(group_hmms2, 'hemopt');
      
    err = totalerror(tmpgroup_matlab, tmpgroup, 4, 0)
    if err>1e-10
      err
      error('unit test for test_hem_hmm_bwd_fwd_mex failed');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% functional test (full)
    clear
    load test_hem_hmm_bwd_fwd_mex_hmms_full.mat
  
    hemopt.sortclusters = 'f';
    hemopt.trials = 5;
    hemopt.emit.covar_type = 'full';
    faceimg = '../demo/face.jpg';
    hemopt
    
    % cluster
    rand('state', 101);
    randn('state', 101);
    
    tic
    [xgroup_hmms2] = vhem_cluster(hmms, 2, 3, hemopt);
    vhem_plot(xgroup_hmms2, faceimg);
    mex_time = toc
  
    rand('state', 101);
    randn('state', 101);
    tic
    hemopt.useMEX = 0
    [group_hmms2] = vhem_cluster(hmms, 2, 3, hemopt);
    vhem_plot(group_hmms2, faceimg);
    matlab_time = toc
  
    tmpgroup        = rmfield(xgroup_hmms2, 'hemopt');
    tmpgroup_matlab = rmfield(group_hmms2, 'hemopt');
      
    err = totalerror(tmpgroup_matlab, tmpgroup, 4, 0)
    if err>1e-10
      err
      error('unit test for test_hem_hmm_bwd_fwd_mex failed');
    end
    
    
  case 3
    %% timing test (diagonal)
    load test_hem_hmm_bwd_fwd_mex_h3ms.mat
    NTrials = 1000;
    
    %% MEX function
    tic
    for n=1:NTrials
      [xL_elbo, xsum_nu_1, xupdate_emit_pr, xupdate_emit_mu, xupdate_emit_Mu, xsum_xi] = ...
        hem_hmm_bwd_fwd_mex(h3m_b.hmm,h3m_r.hmm,T,smooth, 5, 5);
    end
    mex_time = toc
    
    %% MATLAB function
    tic
    for n=1:NTrials
      % storage
      L_elbo           = zeros(Kb,Kr);
      sum_nu_1         = cell(Kb,Kr);
      update_emit_pr   = cell(Kb,Kr);
      update_emit_mu   = cell(Kb,Kr);
      update_emit_Mu   = cell(Kb,Kr);
      sum_xi           = cell(Kb,Kr);
      
      % loop reduced mixture components
      for j = 1 : Kr
        % get this HMM in the reduced mixture
        hmm_r = h3m_r.hmm{j};
        % loop base mixture components
        for i = 1 : Kb
          % get this HMM in the base mixture
          hmm_b = h3m_b.hmm{i};
          
          [L_elbo(i,j) ...
            sum_nu_1{i,j} ...
            update_emit_pr{i,j} ...
            update_emit_mu{i,j} ...
            update_emit_Mu{i,j}  ...
            sum_xi{i,j} ...
            ] = hem_hmm_bwd_fwd(hmm_b,hmm_r,T,smooth, h3m_b_hmm_emit_cache{i});
        end
      end
    end
    matlab_time = toc
    
    % consistency check
    %err = [];
    %err(end+1) = totalerror(xL_elbo, L_elbo, 4, 0);
    %err(end+1) = totalerror(xsum_nu_1, sum_nu_1, 4, 0);
    %err(end+1) = totalerror(xsum_xi, sum_xi, 4, 0);
    %err(end+1) = totalerror(xupdate_emit_pr, update_emit_pr, 4, 0);
    %err(end+1) = totalerror(xupdate_emit_mu, update_emit_mu, 4, 0);
    %err(end+1) = totalerror(xupdate_emit_Mu, update_emit_Mu, 4, 0);
    %if (any(err>1e-10))
    %  err
    %  error('unit test for test_hem_hmm_bwd_fwd_mex failed');
    %end
  
    %% timing test (full) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear
    load test_hem_hmm_bwd_fwd_mex_h3ms_full.mat
    NTrials = 1000;
    
    %% MEX function
    tic
    for n=1:NTrials
      [xL_elbo, xsum_nu_1, xupdate_emit_pr, xupdate_emit_mu, xupdate_emit_Mu, xsum_xi] = ...
        hem_hmm_bwd_fwd_mex(h3m_b.hmm,h3m_r.hmm,T,smooth, maxN, maxN2, logdetCovR, invCovR);
    end
    mex_time = toc
    
    %% MATLAB function
    tic
    for n=1:NTrials
      % storage
      L_elbo           = zeros(Kb,Kr);
      sum_nu_1         = cell(Kb,Kr);
      update_emit_pr   = cell(Kb,Kr);
      update_emit_mu   = cell(Kb,Kr);
      update_emit_Mu   = cell(Kb,Kr);
      sum_xi           = cell(Kb,Kr);
      
      % loop reduced mixture components
      for j = 1 : Kr
        % get this HMM in the reduced mixture
        hmm_r = h3m_r.hmm{j};
        % loop base mixture components
        for i = 1 : Kb
          % get this HMM in the base mixture
          hmm_b = h3m_b.hmm{i};
          
          [L_elbo(i,j) ...
            sum_nu_1{i,j} ...
            update_emit_pr{i,j} ...
            update_emit_mu{i,j} ...
            update_emit_Mu{i,j}  ...
            sum_xi{i,j} ...
            ] = hem_hmm_bwd_fwd(hmm_b,hmm_r,T,smooth, h3m_b_hmm_emit_cache{i});
        end
      end
    end
    matlab_time = toc
end

