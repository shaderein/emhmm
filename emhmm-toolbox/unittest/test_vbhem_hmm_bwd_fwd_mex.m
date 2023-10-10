function test_vbhem_hmm_bwd_fwd_mex(mode)
% test_vbhem_hmm_bwd_fwd_mex - unit test for vbhem_hmm_bwd_fwd_mex
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

% v0.80 - initial version

if (nargin==0)
  mode = 1;
end

switch(mode)
  case 1
    
    %% unit test on MEX function (full covariance)
    % no need to test diag covariance
    clear
    load test_vbhem_hmm_bwd_fwd_mex_h3ms_full.mat
    
    Kb = h3m_b.K;
    Kr = h3m_r.K;
    dim = h3m_b.hmm{1}.emit{1}.nin;
    Sr = length(h3m_r.hmm{1}.emit);
            
    %% MEX function
    [L_elbo, nu_1, update_emit_pr, update_emit_mu, update_emit_M, sum_xi] = ...
         vbhem_hmm_bwd_fwd_mex(h3m_b.hmm,h3m_r.hmm,T,maxN, maxN2, logdetCovPlusDdivlamR, invCovR);
    
    %% MATLAB function
    % storage
    xL_elbo  = zeros(Kb,Kr);
    xnu_1             = cell(Kb,Kr);
    xupdate_emit_pr   = cell(Kb,Kr);
    xupdate_emit_mu   = cell(Kb,Kr);
    xupdate_emit_M    = cell(Kb,Kr);
    xsum_xi           = cell(Kb,Kr);
    
    % cache values
    h3m_b_hmm_emit_cache = make_h3m_b_cache(h3m_b);
      
    err = [];

    % cache the variational parameters of h3m_r
    %% compute phi(rho,beta)
    % loop reduced mixture components
    for j = 1 : Kr
      hmm_r   = h3m_r.hmm{j};
      m       = zeros(dim,Sr);
      W       = zeros(dim,dim,Sr);
      v       = zeros(Sr,1);
      lambda  = zeros(Sr,1);
      logLambdaTilde = zeros(Sr,1);
      
      for k = 1:Sr
        m(:,k) = hmm_r.emit{1, k}.m';
        v(k) = hmm_r.emit{1, k}.v;
        W(:,:,k) = hmm_r.emit{1, k}.W;
        
        lambda(k) = hmm_r.emit{1, k}.lambda;
        logLambdaTilde(k) = h3m_r.hmm{j}.emit{k}.logLambdaTilde;
      end
      
      logATilde  = hmm_r.logATilde;
      logPiTilde = hmm_r.logPiTilde;
      
      % loop base mixture components
      for i = 1 : Kb
        
        hmm_b = h3m_b.hmm{i};
        [xL_elbo(i,j), ...
          xnu_1{i,j}, ...
          xupdate_emit_pr{i,j}, ...
          xupdate_emit_mu{i,j}, ...
          xupdate_emit_M{i,j}, ...
          xsum_xi{i,j} ...
          ] = vbhem_hmm_bwd_fwd_fast(hmm_b,T,m,W,v,lambda,...
          logLambdaTilde,logATilde,logPiTilde, h3m_b_hmm_emit_cache{i});
        
        % calculate errors
        err(end+1) = abs(xL_elbo(i,j)-L_elbo(i,j));
        err(end+1) = sum(abs(xnu_1{i,j}(:) - nu_1{i,j}(:)));
        err(end+1) = sum(abs(xupdate_emit_pr{i,j}(:) - update_emit_pr{i,j}(:)));
        err(end+1) = sum(abs(xupdate_emit_mu{i,j}(:) - update_emit_mu{i,j}(:)));
        err(end+1) = sum(abs(xupdate_emit_M{i,j}(:)  - update_emit_M{i,j}(:)));
        err(end+1) = sum(abs(xsum_xi{i,j}(:) - sum_xi{i,j}(:)));                
      end % end loop base
    end  % end loop reduced                        
    
    % consistency check
    if (any(err>1e-6))
      err      
      error('unit test for test_vhem_hmm_bwd_fwd_mex failed');
    end
    
    
  case 2
    %% funcitonal test (not implemented)
    error('not implemented');
    
  case 3    
    %% timing test (full) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % no need to test diagonal
    clear
    load test_vbhem_hmm_bwd_fwd_mex_h3ms_full.mat
    NTrials = 1000;
    
    Kb = h3m_b.K;
    Kr = h3m_r.K;
    dim = h3m_b.hmm{1}.emit{1}.nin;
    Sr = length(h3m_r.hmm{1}.emit);
            
    %% MEX function
    tic
    for n=1:NTrials
      [L_elbo, nu_1, update_emit_pr, update_emit_mu, update_emit_M, sum_xi] = ...
        vbhem_hmm_bwd_fwd_mex(h3m_b.hmm,h3m_r.hmm,T,maxN, maxN2, logdetCovPlusDdivlamR, invCovR);
    end
    mex_time = toc

    %% MATLAB function
    tic
    for n=1:NTrials
      % storage
      xL_elbo  = zeros(Kb,Kr);
      xnu_1             = cell(Kb,Kr);
      xupdate_emit_pr   = cell(Kb,Kr);
      xupdate_emit_mu   = cell(Kb,Kr);
      xupdate_emit_M    = cell(Kb,Kr);
      xsum_xi           = cell(Kb,Kr);
      
      % cache values
      h3m_b_hmm_emit_cache = make_h3m_b_cache(h3m_b);
      
      % cache the variational parameters of h3m_r
      %% compute phi(rho,beta)
      % loop reduced mixture components
      for j = 1 : Kr
        hmm_r   = h3m_r.hmm{j};
        m       = zeros(dim,Sr);
        W       = zeros(dim,dim,Sr);
        v       = zeros(Sr,1);
        lambda  = zeros(Sr,1);
        logLambdaTilde = zeros(Sr,1);
        
        for k = 1:Sr
          m(:,k) = hmm_r.emit{1, k}.m';
          v(k) = hmm_r.emit{1, k}.v;
          W(:,:,k) = hmm_r.emit{1, k}.W;
          
          lambda(k) = hmm_r.emit{1, k}.lambda;
          logLambdaTilde(k) = h3m_r.hmm{j}.emit{k}.logLambdaTilde;
        end
        
        logATilde  = hmm_r.logATilde;
        logPiTilde = hmm_r.logPiTilde;
        
        % loop base mixture components
        for i = 1 : Kb
          
          hmm_b = h3m_b.hmm{i};
          [xL_elbo(i,j), ...
            xnu_1{i,j}, ...
            xupdate_emit_pr{i,j}, ...
            xupdate_emit_mu{i,j}, ...
            xupdate_emit_M{i,j}, ...
            xsum_xi{i,j} ...
            ] = vbhem_hmm_bwd_fwd_fast(hmm_b,T,m,W,v,lambda,...
            logLambdaTilde,logATilde,logPiTilde, h3m_b_hmm_emit_cache{i});
          
          % calculate errors
          %err(end+1) = abs(xL_elbo(i,j)-L_elbo(i,j));
          %err(end+1) = sum(abs(xnu_1{i,j}(:) - nu_1{i,j}(:)));
          %err(end+1) = sum(abs(xupdate_emit_pr{i,j}(:) - update_emit_pr{i,j}(:)));
          %err(end+1) = sum(abs(xupdate_emit_mu{i,j}(:) - update_emit_mu{i,j}(:)));
          %err(end+1) = sum(abs(xupdate_emit_M{i,j}(:)  - update_emit_M{i,j}(:)));
          %err(end+1) = sum(abs(xsum_xi{i,j}(:) - sum_xi{i,j}(:)));
        end % end loop base
      end  % end loop reduced
    end
    
    matlab_time = toc
end

