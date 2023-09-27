function h3m_rs2dim = hem_coh3m_c_step(h3m_rs,h3m_bs,mopt,emopt) 
% modified function for co-clustering
%
% h3m_r = initialization for the reduced mixture output
% h3m_b = base mixture input
% mopt  = options 
%       .K     = number of mixtures in the reduced model
%       .Nv     = number of virtual samples (this is multiplied by Kb)
%       .tau   = length of virtual sequences
%       .termmode  = how to terminate EM
%       .termvalue = when to terminate the EM
% emopt = info about number of stimuli & subject
% h3m.K
% h3m.hmm = {hmm1 hmm2 ... hmmK}
% h3m.omega = [omega1 omega2 ... omegaK]
% hmm1.A hmm1.emit hmm1.prior
%
% ---
% H3M Toolbox 
% Copyright 2012, Emanuele Coviello, CALab, UCSD

% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2019-07-13
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

% 2016-12-09: ABC - added output for emission population size
% 2017-04-25: ABC - small bug fix when generate component is handled.
% 2017-06-06: ABC - put M-step into a separate function
% 2017-08-07: ABC - caching to make E-step faster
% 2017-08-08: ABC - added MEX implementation for E-step, and useMEX option
% 2017-08-17: ABC - added "full" covariance matrices

% 2020-03-12: v0.77 - initial version
% 2020-08-06: v0.80 - fix bug w/ empty HMMs

tic;

h3m_rs2dim.h3m_rs = h3m_rs;
h3m_rs2dim.LogL = [];
h3m_rs2dim.K = h3m_rs(1).K;

num_iter = 0;

% number of components in the base and reduced mixtures
Kb = emopt.subjectsize;
Kr = h3m_rs2dim.K;

% number of states 
Ns = [];
Nrs = [];
for i=1:emopt.stimulisize
    N = size(h3m_bs(i).hmm{1}.A,1);
    Nr = size(h3m_rs(i).hmm{1}.A,1);
    Ns = [Ns; N];
    Nrs = [Nrs; Nr];
end

% number of mixture components in each state emission probability
Ms = [];
for i=1:emopt.stimulisize
    M = h3m_bs(i).hmm{1}.emit{1}.ncentres;
    Ms = [Ms; M];
end

% length of the virtual sample sequences
T = mopt.tau;

% dimension of the emission variable
dims = [];
for i=1:emopt.stimulisize
    dim = h3m_bs(i).hmm{1}.emit{1}.nin;
    dims = [dims; dim];
end

% number of virtual samples 
virtualSamples = mopt.Nv;

% correct
% ABC 2020-08-06: BUG FIX v0.80: compute N_i for each stimuli
% it can be different because some trials may be missing for different stimuli
% and omega=0 for those trials
%N_i = virtualSamples * h3m_bs(1).omega * Kb;
Nk_i = {};
for k=1:emopt.stimulisize
  tmp = virtualSamples * h3m_bs(k).omega * Kb;
  Nk_i{k} = tmp(:);
end

% regularization value
if isfield(mopt,'reg_cov')
    reg_cov = mopt.reg_cov;
else
    reg_cov = 0;
end

% minimum number of iterations
if ~isfield(mopt,'min_iter')
    mopt.min_iter = 0;
end

if ~isfield(mopt,'story')
    mopt.story = 0;
end

% timing
here_time = tic;
if isfield(mopt,'start_time')
    a_time = mopt.start_time;
else
    a_time = here_time;
end

% save the H3Ms in each iteration
if mopt.story
    story = {};
    time = [];
end
if mopt.story
    story{end+1} = h3m_rs;
    time(end+1) = toc(a_time);
end

% regularize
for j = 1 : Kr
  for i = 1 : length(Nrs)
      for n = 1 : Nrs(i)  % bug fix (was N before)
        switch(h3m_rs(i).hmm{j}.emit{n}.covar_type)
          case 'diag'        
            h3m_rs(i).hmm{j}.emit{n}.covars = h3m_rs(i).hmm{j}.emit{n}.covars + reg_cov;
          case 'full'
            h3m_rs(i).hmm{j}.emit{n}.covars = h3m_rs(i).hmm{j}.emit{n}.covars + reg_cov*eye(dim);
        end
      end
  end
end

switch mopt.inf_norm
  case ''
    inf_norm = 1;
  case 'n'
    inf_norm = virtualSamples / Kb;
  case {'tn' 'nt'}
    inf_norm = T * virtualSamples / Kb;
  case {'t'}
    inf_norm = T;  
end

smooth =  mopt.smooth;

% 2017-08-08 - ABC - handle MEX 
if ~isfield(mopt, 'useMEX')
  mopt.useMEX = 1;
end

canuseMEX = 0;
if (mopt.useMEX)
  if (exist('hem_hmm_bwd_fwd_mex')==3)
    canuseMEX = 1;
    for i=1:length(h3m_bs)
        for n=1:length(h3m_bs(i).hmm)
          for j=1:length(h3m_bs(i).hmm{n}.emit)
            if h3m_bs(i).hmm{n}.emit{j}.ncentres ~= 1
              canuseMEX = 0;
              warning('hem_hmm_bwd_fwd_mex: ncentres should be 1');
            end
            %if ~strcmp(h3m_b.hmm{n}.emit{j}.covar_type, 'diag')
            %  canuseMEX = 0;
            %  warning('hem_hmm_bwd_fwd_mex: only diag is supported');
            %end
            if canuseMEX == 0
              break;
            end
          end
          if canuseMEX == 0
            break
          end
        end
    end
    
  else
    canuseMEX = 0;
    warning('hem_hmm_bwd_fwd_mex: MEX implementation not found');
  end
end

if (canuseMEX)
  % get maxN and maxN2
  maxNs = [];
  maxN2s = [];
  for i=1:emopt.stimulisize
      maxN  = max(cellfun(@(x) length(x.prior), h3m_bs(i).hmm));
      maxN2 = max(cellfun(@(x) length(x.prior), h3m_rs(i).hmm));
      maxNs = [maxNs; maxN];
      maxN2s = [maxN2s; maxN2];
  end
  
else
  % 2017-08-07 - ABC - Caching/reshaping for faster E-step
  h3m_b_hmm_emit_cache = cell(emopt.stimulisize,Kb);
  for k=1:emopt.stimulisize
      for i=1:Kb
        g3m_b = h3m_bs(k).hmm{i}.emit;
        BKall = cellfun(@(x) x.ncentres, g3m_b);
        BKmax = max(BKall);

        % does not support different number of centres (for speed)
        if ~all(BKall==BKmax)
          error('different number of centres is not supported.');
        end

        gmmBcentres_tmp = cellfun(@(x) x.centres, g3m_b, ...
          'UniformOutput', false);
        h3m_b_hmm_emit_cache{k,i}.gmmBcentres = cat(3,gmmBcentres_tmp{:});  % BKmax x dim x N

        % extract all the covars
        gmmBcovars_tmp = cellfun(@(x) x.covars, g3m_b, ...
          'UniformOutput', false);
        switch (g3m_b{1}.covar_type)
          case 'diag'
            h3m_b_hmm_emit_cache{k,i}.gmmBcovars = cat(3,gmmBcovars_tmp{:});  % BKmax x dim x N
          case 'full'
            h3m_b_hmm_emit_cache{k,i}.gmmBcovars = cat(4,gmmBcovars_tmp{:});  % dim x dim x BKmax x N
            % TODO: add cache for logdet(gmmBcovars)
        end

        % extract all priors
        gmmBpriors_tmp = cellfun(@(x) x.priors', g3m_b, ...
          'UniformOutput', false);
        h3m_b_hmm_emit_cache{k,i}.gmmBpriors = cat(2,gmmBpriors_tmp{:});  % BKmax x N
      end
  end
end

% initialize arrays for E-step
  L_elbo           = zeros(emopt.stimulisize,Kb,Kr);
  nu_1             = cell(emopt.stimulisize,Kb,Kr);
  update_emit_pr   = cell(emopt.stimulisize,Kb,Kr);
  update_emit_mu   = cell(emopt.stimulisize,Kb,Kr);
  update_emit_M    = cell(emopt.stimulisize,Kb,Kr);
  sum_xi           = cell(emopt.stimulisize,Kb,Kr);

% start looping variational E step and M step
while 1

  %%%%%%%%%%%%%%%%%%%%
  %%%    E-step    %%%
  %%%%%%%%%%%%%%%%%%%%
    
  if (canuseMEX)
    %% MEX implementation
    
    switch(h3m_bs(1).hmm{1}.emit{1}.covar_type)
      case 'diag'
        for k=1:emopt.stimulisize
            [L_elbo(k,:,:), nu_1(k,:,:), update_emit_pr(k,:,:), update_emit_mu(k,:,:), update_emit_M(k,:,:), sum_xi(k,:,:)] = ...
          hem_hmm_bwd_fwd_mex(h3m_bs(k).hmm,h3m_rs(k).hmm,T,smooth, maxN(k), maxN2(k));
        end
      case 'full'
        for k=1:emopt.stimulisize
            % cache logdets and inverse covariances
            logdetCovR = cell(1,Kr);
            invCovR    = cell(1,Kr);
            for j=1:Kr
               % does not support ncentres>1
              logdetCovR{j} = zeros(1,length(h3m_rs(k).hmm{j}.prior));
              invCovR{j}    = zeros(dim, dim, length(h3m_rs(k).hmm{j}.prior));
              for n=1:length(h3m_rs(k).hmm{j}.prior)
                logdetCovR{j}(n)  = log(det(h3m_rs(k).hmm{j}.emit{n}.covars)); 
                invCovR{j}(:,:,n) = inv(h3m_rs(k).hmm{j}.emit{n}.covars);
              end
            end
            [L_elbo(k,:,:), nu_1(k,:,:), update_emit_pr(k,:,:), update_emit_mu(k,:,:), update_emit_M(k,:,:), sum_xi(k,:,:)] = ...
            hem_hmm_bwd_fwd_mex(h3m_bs(k).hmm,h3m_rs(k).hmm,T,smooth, maxNs(k), maxN2s(k), logdetCovR, invCovR);
        end
    end
    
    % consistency check
    if 0
      err = [];
      err(end+1) = totalerror(xL_elbo, L_elbo);
      err(end+1) = totalerror(xnu_1, nu_1);
      err(end+1) = totalerror(xsum_xi, sum_xi);
      err(end+1) = totalerror(xupdate_emit_pr, update_emit_pr);
      err(end+1) = totalerror(xupdate_emit_mu, update_emit_mu);
      err(end+1) = totalerror(xupdate_emit_M, update_emit_M);
      if any(err>1e-6)
        warning('mismatch');
        keyboard
      end
    end
    
  else
    %% MATLAB implementation
    % iterate through stimuli
    for k = 1 : emopt.stimulisize
        % loop reduced mixture components
        for j = 1 : Kr

          % get this HMM in the reduced mixture
          hmm_r = h3m_rs(k).hmm{j};

          % loop base mixture components
          for i = 1 : Kb

            % get this HMM in the base mixture
            hmm_b = h3m_bs(k).hmm{i};

            % run the Hierarchical bwd and fwd recursions. This will compute:
            % the ELBO on E_M(b)i [log P(y_1:T | M(r)i)]
            % the nu_1^{i,j} (rho)
            % the sum_gamma  [sum_t nu_t^{i,j} (rho,gamma)] sum_m c.. nu..
            % the sum_gamma  [sum_t nu_t^{i,j} (rho,gamma)] sum_m c.. nu.. mu ..
            % the sum_gamma  [sum_t nu_t^{i,j} (rho,gamma)] sum_m c.. nu.. M .. which accounts for the mu mu' and the sigma
            % [the last 3 are not normalized, of course]
            % the sum_t xi_t^{i,j} (rho,sigma)     (t = 2 ... T)
            [L_elbo(k,i,j) ...
              nu_1{k,i,j} ...
              update_emit_pr{k,i,j} ...
              update_emit_mu{k,i,j} ...
              update_emit_M{k,i,j}  ...
              sum_xi{k,i,j} ...
              ] = hem_hmm_bwd_fwd(hmm_b,hmm_r,T,smooth, h3m_b_hmm_emit_cache{k,i});

            % consistency check
            if 0
              [xL_elbo ...
                xnu_1 ...
                xupdate_emit_pr ...
                xupdate_emit_mu ...
                xupdate_emit_M  ...
                xsum_xi ...
                ] = hem_hmm_bwd_fwd_OLD(hmm_b,hmm_r,T,smooth);
              err = [];
              err(end+1) = abs(xL_elbo-L_elbo(i,j));
              err(end+1) = sum(abs(xnu_1(:) - nu_1{i,j}(:)));
              err(end+1) = sum(abs(xupdate_emit_pr(:) - update_emit_pr{i,j}(:)));
              err(end+1) = sum(abs(xupdate_emit_mu(:) - update_emit_mu{i,j}(:)));
              err(end+1) = sum(abs(xupdate_emit_M(:)  - update_emit_M{i,j}(:)));
              err(end+1) = sum(abs(xsum_xi(:) - sum_xi{i,j}(:)));
              if any(err>1e-6)
                warning('mismatch')
                keyboard;
              end

            end
            % end loop base
          end

          % end loop reduced
       end
    end
  end
  
  %% ABC - if first iteration, sort the reduced components to be more consistent
  %{
  if (num_iter==0)
    % L_elbo [stimuli x base x reduced]
    newinds = match_clusters(L_elbo);
    
    % permute
    for k=1:emopt.stimulisize
      knewinds = newinds(k,:);
      L_elbo(k,:,:) = L_elbo(k,:,knewinds);            
      nu_1(k,:,:) = nu_1(k,:,knewinds);
      update_emit_pr(k,:,:) = update_emit_pr(k,:,knewinds);
      update_emit_mu(k,:,:) = update_emit_mu(k,:,knewinds);
      update_emit_M(k,:,:) = update_emit_M(k,:,knewinds);
      sum_xi(k,:,:)        = sum_xi(k,:,knewinds);
      h3m_rs(k).hmm = h3m_rs(k).hmm(knewinds);
    end
  end
  %}
  
    
  %% ABC - checked
  % compute the z_ij
  % this is not normalized ...
  % log_Z = ones(Kb,1) * log(h3m_r.omega) + diag(N_i) * L_elbo;
  L_elbo = L_elbo / (inf_norm);
  log_Z = zeros(Kb,Kr);
  for k=1:emopt.stimulisize
     % Note: the Nk_i term contains 0 for base HMMs with omega=0
     log_Z = log_Z + (Nk_i{k}(:) * ones(1,Kr)) .* reshape(L_elbo(k,:,:),[Kb,Kr]);
  end 
  log_Z = log_Z / emopt.stimulisize;
  log_Z = log_Z + ones(Kb,1) * log(h3m_rs(1).omega);
  
  % normalize Z
  Z = exp(log_Z - logtrick(log_Z')' * ones(1,Kr));
  
  % compute the elbo to total likelihood ...  
  %     ll = Z .* (log_Z - log(Z));
  %
  %     ll(isnan(ll)) = 0;
  %
  %     % new_LogLikelihood = sum(sum(ll));
  %     new_LogLikelihood_foo = ones(1,Kb) * ll * ones(Kr,1);
    
  % these should be the same
  new_LogLikelihood = sum(logtrick(log_Z')');
    
  % update the log likelihood in the reduced mixture
  old_LogLikelihood = h3m_rs(1).LogL;
  for i=1:emopt.stimulisize
      h3m_rs(i).LogL = new_LogLikelihood;
      %h3m_rs(i).LogL(end+1) = new_LogLikelihood;
      h3m_rs(i).Z = Z;
  end
    
  % check whether to continue with a new iteration or to return
  stop = 0;
  if num_iter > 1
    changeLL = (new_LogLikelihood - old_LogLikelihood) / abs(old_LogLikelihood);
  else
    changeLL = inf;
  end
  
  if (mopt.verbose > 1)
    fprintf('iter=%d; LL=%g (pL=%g)\n', num_iter, new_LogLikelihood, changeLL);
  end
  
  if changeLL<0
    %warning('The change in log likelihood is negative!!!')
    fprintf('[N]');
  end
  
  % change in LL is below threshold?
  switch mopt.termmode
    case 'L'
      if (changeLL < mopt.termvalue)
        stop = 1;
      end
  end
  
  % max iteration is reached?
  if (num_iter > mopt.max_iter)
    stop = 1;
  end
        
  if isnan(new_LogLikelihood)
    fprintf('[NaN]');
    stop = 1;
  end
  
  % stop and reached minimum number of iterations?
  if (stop) && (num_iter >= mopt.min_iter)
    % then exit the EM loop
    break
  end
  
  num_iter = num_iter + 1;  
  old_LogLikelihood = new_LogLikelihood;
    
  %%%%%%%%%%%%%%%%%%%%
  %%%    M-step    %%%
  %%%%%%%%%%%%%%%%%%%%
  % compute new parameters

  % make a copy first
  h3m_rs_new = h3m_rs;
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%  re-estimation of the component weights (omega)
  omega_new = (ones(1,Kb) / Kb) * Z;
  h3m_r_new.omega = omega_new;
  
  for k = 1:emopt.stimulisize
      %% for each reduced component j
      for j = 1:Kr    
        % get statistics for j-th component
        estats.Zomega  = Z(:,j) .* Nk_i{k}(:);
        estats.nu_1    = nu_1(k,:,j);
        estats.sum_xi  = sum_xi(k,:,j);
        estats.emit_pr = update_emit_pr(k,:,j);
        estats.emit_mu = update_emit_mu(k,:,j);
        estats.emit_M  = update_emit_M(k,:,j);

        % setup constants
        estats.M   = M;
        estats.N2  = size(h3m_rs(k).hmm{j}.A,1);
        estats.dim = dim;
        estats.Kb  = Kb;
        estats.reg_cov = reg_cov;

        % run M-step on the j-th component
        h3m_rs_new(k).hmm{j} = hem_mstep_component(h3m_rs_new(k).hmm{j}, estats, mopt);
      end
  end
  
  %% check for degenerate components %%%%%%%%%%%%%%%%%
        
  % ABC: 2020-08-06: use a threshold
  OMEGA_THRESH = 0;
  %OMEGA_THRESH = 1e-5;
  % if a component gets 0 prior ...
  ind_zero = find(h3m_r_new.omega <= OMEGA_THRESH);
  for k=1:emopt.stimulisize
      for i_z = ind_zero
        fprintf('[P]');
        %% ABC BUG FIX: Nr -> Nrs(k), N -> Ns(k)        
        h3m_rs_new(k) = hem_fix_degenerate_component(h3m_rs_new(k), i_z, Ns(k), Nrs(k));
      end
  end
  
  % ABC: 2018-02-13
  % check if an HMM state has zero counts
  % ABC: 2020-08-06: use a threshold 
  VCOUNT_THRESH = 0;
  %VCOUNT_THRESH = 1e-5;
  for k = 1 : emopt.stimulisize
      for j = 1 : Kr
        ind_zero = find(h3m_rs_new(k).hmm{j}.stats.emit_vcounts <= VCOUNT_THRESH);
        if ~isempty(ind_zero)
          for i_z = ind_zero            
            fprintf('[S]');
            h3m_rs_new(k).hmm{j} = hem_fix_degenerate_hmm(h3m_rs_new(k).hmm{j}, i_z);
          end
        end
      end
  end
    
  % check if an emission has zero prior (only applicable for GMM emission with >1 components
  if (h3m_rs_new(1).hmm{1}.emit{1}.ncentres > 1)
    for k = 1 : emopt.stimulisize
        for j = 1 : Kr
          N2 = size(h3m_rs(k).hmm{j}.A,1);      
          for n = 1 : N2
            % if some of the emission prob is zero, replace it ...
            ind_zero = find(h3m_rs_new(k).hmm{j}.emit{n}.priors == 0);
            for i_z = ind_zero              
              %fprintf('[E]');
              h3m_rs_new(k).hmm{j}.emit{n} = hem_fix_degenerate_emission(h3m_rs_new(k).hmm{j}.emit{n}, i_z);
            end
          end
        end
    end
  end
  
        
  %% update the model %%%%%%%%%%%%%%%%%%%%
  h3m_rs = h3m_rs_new;
  h3m_rs2dim.LogL = h3m_rs(1).LogL;
  h3m_rs2dim.h3m_rs = h3m_rs;
  
  % add to story
  if mopt.story
    story{end+1} = h3m_rs;
    time(end+1) = toc(a_time);
  end  
end

for k=1:emopt.stimulisize
    h3m_rs2dim.elapsed_time = toc(here_time);
    if mopt.story
      h3m_rs2dim.story = story;
      h3m_rs2dim.time = time;
    end
end

% end of the function
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s] = logtrick(lA)
% logtrick - "log sum trick" - calculate log(sum(A)) using only log(A) 
%
%   s = logtrick(lA)
%
%   lA = column vector of log values
%
%   if lA is a matrix, then the log sum is calculated over each column
% 

[mv, mi] = max(lA, [], 1);
temp = lA - repmat(mv, size(lA,1), 1);
cterm = sum(exp(temp),1);
s = mv + log(cterm);
end



