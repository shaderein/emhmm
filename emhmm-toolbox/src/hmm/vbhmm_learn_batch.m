function [hmm_batch, L_batch, vbopt2] = vbhmm_learn_batch(data_batch, K, vbopt)
% vbhmm_learn_batch - Learn HMMs on a batch of subjects
%
%   [hmm_batch, L_batch, vbopt2] = vbhmm_learn_batch(data_batch, K, vbopt)
%
%  INPUTS
%    data_batch = Sx1 cell array, one cell for each subject.
%                 Each cell entry contains a cell array of N fixation sequences.
%                 i.e., the i-th trial of the s-th subject is:
%                       data{s}{i} = [TxD] matrix, where each row is a fixation (x,y) or (x,y,t)
%                         N = number of sequences, T = sequence length.
%                         D = dimension: 2 for fixation location (x,y); 3 for location & duration (x,y,d). Duration should be in milliseconds.
%                 Note: the sequence length T, and number of sequences N can be different for each subject.
%
%      K = scalar: number of hidden states (clusters)
%          vector: automatically selects the number of hidden states from the given values. 
%                  The model K with highest log-likelihood is selected.
%
%  vbopt =  VB options structure. See vbhmm_learn for more details.
%
%  OUTPUTS
%    hmm_batch = 1xS cell array containing an HMM for each subject
%      L_batch = 1xS matrix with log-likelihoods for each subject
%       vbopt2 = the options used if estimating hyperparameters automatically
%
% NOTES
%   To improve the speed, this function uses parallel for loops (parfor)
%   from the Parallel Computing Toolbox.
%  
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-08-04
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

% VERSIONS
% 2018-11-23: v0.74 - skip if empty data
% 2019-01-12: v0.75 - use parfor to improve speed
% 2019-02-20: v0.75 - make verbose output better
% 2019-02-21: v0.75 - parfor for batch hyp learning
% 2021-07-22: v0.78 - support for fixed_rois with learn_hyps_batch

% v0.75 - verbose check
if ~isfield(vbopt, 'verbose')
  vbopt.verbose = 1;
end


%% no hyperparameter estimation
if ~isfield(vbopt, 'learn_hyps_batch') || (~iscell(vbopt.learn_hyps_batch) && vbopt.learn_hyps_batch==0)
  
  N = length(data_batch);
  hmm_batch = cell(1,N);
  L_batch = zeros(1,N);
  
  % loop over each subject, and estimate the HMM
  parfor i=1:N
    
    fprintf('=== running Subject %d ===\n', i);
    
    % 2018-11-23: v0.74 - skip if empty data
    if isempty(data_batch{i})
      fprintf('** no data...skipping **\n');
      hmm_batch{i} = [];
      L_batch(i) = -inf;
      
    else
      % v0.75 - tag each verbose output
      vbopt2 = vbopt;
      vbopt2.verbose_prefix = sprintf('[Subject %d] ', i);
      
      hmm_batch{i} = vbhmm_learn(data_batch{i}, K, vbopt2);
      L_batch(i) = hmm_batch{i}.LL;
      if vbopt.showplot
        set(gcf, 'name', sprintf('Subject %d', i));
        drawnow
      end
    end
  end
    
  vbopt2 = vbopt;

else 
  
  %% do hyperparameter estimation (in batch)
  vbopt = setdefault(vbopt, 'minimizer', 'minimize-bfgs');
  
  %% set default hyps
  if ~iscell(vbopt.learn_hyps_batch)
    if (vbopt.learn_hyps_batch == 1)
      vbopt.learn_hyps_batch = {'alpha0', 'epsilon0', 'v0', 'beta0', 'W0', 'mu0'};
      
      % v0.78: fixed ROIs - auto learning only alpha0 and epsilon0
      if isfield(vbopt, 'fixed_rois') && iscell(vbopt.fixed_rois)
        vbopt.learn_hyps_batch = {'alpha0', 'epsilon0'};
      end
      
    else
      error('learn_hyps_batch not set properly');
    end
  else
    % v0.78: if using fixed ROIs, check if valid hyps set
    if isfield(vbopt, 'fixed_rois') && iscell(vbopt.fixed_rois)
      for i=1:length(vbopt.learn_hyps_batch)
        if ~any(strcmp(vbopt.learn_hyps_batch{i}, {'alpha0', 'epsilon0'}))
          error('fixed ROIs: cannot set hyp %s as learnable', vbopt.learn_hyps_batch{i})
        end
      end   
    end
  end
  learn_hyps_batch = vbopt.learn_hyps_batch;
  
  % remove option for individual learning
  if isfield(vbopt, 'learn_hyps')
    warning('vbopt.learn_hyps_batch and vbopt.learn_hyps were both set. Using only vbopt.learn_hyps_batch.');
    vbopt = rmfield(vbopt, 'learn_hyps');
  end
    
  %% get the hyps for optimization
  hypinfo = get_hypinfo(learn_hyps_batch, vbopt);
  
  %% get the initial HMMs
  inithmms = cell(1,length(data_batch));
  % v0.75 - speedup with parfor
  parfor n=1:length(data_batch)    
    vbopt2 = vbopt;
    vbopt2.keep_suboptimal_hmms = 1;
    vbopt2.showplot = 0;  % disable plotting for initial HMMs
    vbopt2.verbose_prefix = sprintf('[Subject %d] ', n);
    
    fprintf('=== learning initial HMM for subject %d ===\n', n);
    inithmms{n} = vbhmm_learn(data_batch{n}, K, vbopt2);
  end
 
  
  fprintf('=== Optimizing hyperparameters ===\n');
  
  %% initial parameter vector, and the objective function
  initX = init_hyp(vbopt, hypinfo);
  myf   = @(X) vbhmm_grad_batch_parfor(data_batch, inithmms, X, hypinfo);
  %myf   = @(X) vbhmm_grad_batch(data_batch, inithmms, X, hypinfo);
  
  switch(vbopt.minimizer)
    case 'fminunc'
      %% do optimization using fminunc      
      % set optimization options
      options = optimoptions(...
        'fminunc', ...
        'Algorithm', 'trust-region', ...
        'GradObj',   'on', ...
        'Display', 'off', ...
        'TolFun',    inithmms{1}.vbopt.minDiff, ...
        'OutputFcn', @outfun);
      %    'Display',   'Iter-Detailed',
      
      % do the optimization
      [opt_transhyp, opt_nL, exitflag] = fminunc(myf, initX, options);
      
      %opt_transhyp'
      %opt_nL
      
      % check for problem in exitflag
      if (exitflag==0)
        warning('not enough iterations to fully optimize');
      elseif (exitflag<0)
        error('error optimizing the function');
      end
      
    case {'minimize-lbfgs', 'minimize-bfgs', 'minimize-cg'}
      %% do the optimization using minimize_new (LBFGS)
      p.length = 100;
      if strcmp(vbopt.minimizer, 'minimize-lbfgs')
        p.method = 'LBFGS';
      elseif strcmp(vbopt.minimizer, 'minimize-cg')
        p.method = 'CG';
      else
        p.method = 'BFGS';
      end
      p.verbosity = 2; %inithmms{1}.vbopt.verbose;
      [opt_transhyp, opt_Ls] = minimize_new(initX, myf, p);
      opt_nL = opt_Ls(end);
      fprintf('final LL=%g\n', opt_nL);
      
      %[opt_transhyp, opt_Ls] = minimize(initX, myf, 100);
    
  otherwise
    error('bad minimizer specified');
  end
  
  %% run EM to get the final model again
  % (since fminunc doesn't return results from function evaluations)
  fprintf('final run: ');
  [finalnL, finalddL, tmp, hmm_batch] = vbhmm_grad_batch_parfor(data_batch, inithmms, opt_transhyp, hypinfo);
  %[finalnL, finalddL, tmp, hmm_batch] = vbhmm_grad_batch(data_batch, inithmms, opt_transhyp, hypinfo);
  
  L_batch = zeros(1,length(data_batch));
  for n=1:length(data_batch)
    L_batch(n) = hmm_batch{n}.LL;
    
    %% add info about learning hyps
    %hmm_batch{n}.learn_hyps.hypinfo      = hypinfo;
    %hmm_batch{n}.learn_hyps.opt_transhyp = opt_transhyp;
    %hmm_batch{n}.learn_hyps.opt_L        = -opt_L;  % opt_L is actually the negative Log-likelihood
    %hmm_batch{n}.learn_hyps.inithmm      = vbopt2.inithmm;
  end
  
  %% final processing
  vbopt2 = hmm_batch{1}.vbopt;
  
  % dump optimized hyps
  fprintf('\nOptimized hyperparameters: {');
  for k=1:length(hypinfo)
    myname = hypinfo(k).optname;
    fprintf('%s=%s; ', myname, mat2str(vbopt2.(myname), 4));
  end
  fprintf('}\n');
  
  % post processing
  % - reorder based on cluster size
  if ~isempty(vbopt2.sortclusters)
    hmm_old_batch = hmm_batch;
    for i=1:length(data_batch)
      hmm_batch{i} = vbhmm_standardize(hmm_batch{i}, vbopt2.sortclusters);
    end
  end
  
  % show plots
  if (vbopt.showplot)
    for i=1:length(data_batch)
      vbhmm_plot(hmm_batch{i}, data_batch{i}, vbopt.bgimage);
      set(gcf, 'name', sprintf('Subject %d', i));    
      drawnow
    end
  end
  
  % show subject
  if 0
    vbhmm_plot(opt_hmm, mydata, faceimg);
    figure,
    subplot(2,1,1)
    vbhmm_plot_compact(hmm, faceimg);
    title(sprintf('original LL=%g', hmm.LL));
    subplot(2,1,2)
    vbhmm_plot_compact(opt_hmm, faceimg);
    title(sprintf('hyp optimized LL=%g', opt_hmm.LL));
    
    %% plot new m
    if any(strcmp(optname, 'mu'))
      figure
      img = imread(faceimg);
      plot_emissions(mydata, opt_hmm.gamma, opt_hmm.pdf, img)
      hold on
      legs(1) = plot(vbopt3.mu(1), vbopt3.mu(2), 'kx', 'Markersize', 10);
      legs(2) = plot(vbopt.mu(1), vbopt3.mu(2), 'ko', 'Markersize', 10);
      meanfix = mean(cat(1,mydata{:}));
      legs(3) = plot(meanfix(1), meanfix(2), 'k+', 'Markersize', 10);
      hold off
      legend(legs, {'optimized', 'old', 'mean'});
    end
  end

  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% function to modify vbopt using parameter vector X %%%%%%%%%%%%%%%%%%%%%%
function vbopt2 = set_vbopt(X, vbopt, hypinfo)
vbopt2 = vbopt;
ind = 1;
for i=1:length(hypinfo)
  myinds = ind:(ind+hypinfo(i).hypdims-1);
  transhyp = hypinfo(i).hyptrans(X(myinds));  % transform from opt-space to hyp-space
  vbopt2.(hypinfo(i).optname) = transhyp;
  ind = ind+hypinfo(i).hypdims;
end


%% function for calcluating the objective and gradient %%%%%%%%%%%%%%%%%%%
function [nL, dnL, H, hmm_batch] = vbhmm_grad_batch(data_batch, inithmm_batch, X, hypinfo)
% X is the parameter vector

nL = 0;
dnL = [];
H = [];
hmm_batch = cell(1,length(data_batch));

for n=1:length(data_batch)
  % get the used vbopt, and set new parameters
  inithmm = inithmm_batch{n};
  vbopt2  = inithmm.vbopt;
  vbopt2  = set_vbopt(X, vbopt2, hypinfo);
  
  %fprintf('%d', n);
  fprintf('.');
  
  Q = length(inithmm.suboptimal_hmms);
  qhmms = cell(1,Q);
  qnLs  = zeros(1,Q);
  for q=1:Q
%    fprintf('.');
    % run EM using inithmm batches
    vbopt2.initmode = 'inithmm';
    vbopt2.inithmm  = inithmm.suboptimal_hmms{q};
    K = length(vbopt2.inithmm.pdf);
    if (nargout>1)
      vbopt2.calc_LLderiv = 1;
    end
    
    qhmms{q} = vbhmm_em(data_batch{n}, K, vbopt2);
    
    % since K can be different for each init,
    % correct for multiple equivalent parameterizations
    qnLs(q) = -qhmms{q}.LL - gammaln(K+1);      
  end
  
  % pick the best (minimize)
  [bestL, bestind] = min(qnLs);
  %fprintf('(%d)', bestind);
  %fprintf(' ');
  
  
  mynL = qnLs(bestind);
  hmm = qhmms{bestind};
  
  % update LL
  nL = nL + mynL;
  
  if (nargout>1)
    mydL = [];
    % append gradients (negate since we actually want to maximize)
    for i=1:length(hypinfo)
      mydL = [mydL; -hmm.dLL.(hypinfo(i).derivname)(:)];
    end
    
    % accumulate gradient
    if isempty(dnL)
      dnL = mydL;
    else
      dnL = dnL + mydL;
    end
  end
  
  if (nargout>3)
    hmm_batch{n} = hmm;
    hmm_batch{n}.vbopt = vbopt2;
  end
end


%fprintf('\n');
%fprintf('.');
fprintf('|');

% normalize by batch length (for stopping criteria)
nL  = nL / length(data_batch);
dnL = dnL / length(data_batch);



%% function for calcluating the objective and gradient %%%%%%%%%%%%%%%%%%%
% v0.75 - use parfor for gradient calculations
function [nL, dnL, H, hmm_batch] = vbhmm_grad_batch_parfor(data_batch, inithmm_batch, X, hypinfo)
% X is the parameter vector

nL = 0;
dnL = [];
H = [];
N = length(data_batch);
hmm_batch = cell(1,length(N));

if (nargout>1)
  do_deriv = 1;
else
  do_deriv = 0;
end

% build combinations of (n,q)
Qs = zeros(1,N);
nqlen = 0;
for n=1:N
  Qs(n) = length(inithmm_batch{n}.suboptimal_hmms);
end
nqlen = sum(Qs);
nqs = zeros(nqlen,2);
ind = 1;
for n=1:N
  nqs(ind+(0:Qs(n)-1),:) = [n*ones(Qs(n),1), (1:Qs(n))'];
  ind = ind+Qs(n);
end

% storage
all_qhmms = cell(1,nqlen);
all_qnLs  = zeros(1,nqlen);
all_vbopt2 = cell(1,nqlen);

% for each (n,q) combination
parfor nq=1:nqlen
  % get (n,q)
  n = nqs(nq,1);
  q = nqs(nq,2);
  
  %fprintf('.');
  
  % get the used vbopt, and set new parameters
  inithmm = inithmm_batch{n};
  vbopt2  = inithmm.vbopt;
  vbopt2  = set_vbopt(X, vbopt2, hypinfo);
  
  % set the inithmm
  vbopt2.initmode = 'inithmm';
  vbopt2.inithmm  = inithmm.suboptimal_hmms{q};
  K = length(vbopt2.inithmm.pdf);
  if (do_deriv)
    vbopt2.calc_LLderiv = 1;
  end
    
  % run EM on (n,q)
  qhmm = vbhmm_em(data_batch{n}, K, vbopt2);
    
  % since K can be different for each init,
  % correct for multiple equivalent parameterizations
  qnL = -qhmm.LL - gammaln(K+1);
  
  % store
  all_qhmms{nq} = qhmm;
  all_qnLs(nq)  = qnL;
  all_vbopt2{nq} = vbopt2;
end

% aggregate the results
ind_start = 0;
for n=1:N
  myinds = ind_start + (1:Qs(n));
  
  % pick the best (minimize)
  [bestL, ibestind] = min(all_qnLs(myinds));
  bestind = myinds(ibestind);  
  mynL = all_qnLs(bestind);
  hmm = all_qhmms{bestind};
  vbopt2 = all_vbopt2{bestind};
  
  % update LL
  nL = nL + mynL;
  
  if (nargout>1)
    mydL = [];
    % append gradients (negate since we actually want to maximize)
    for i=1:length(hypinfo)
      mydL = [mydL; -hmm.dLL.(hypinfo(i).derivname)(:)];
    end
    
    % accumulate gradient
    if isempty(dnL)
      dnL = mydL;
    else
      dnL = dnL + mydL;
    end
  end
  
  if (nargout>3)
    hmm_batch{n} = hmm;
    hmm_batch{n}.vbopt = vbopt2;
  end
  
  ind_start = ind_start + Qs(n);
end

fprintf('|');

% normalize by batch length (for stopping criteria)
nL  = nL / N;
dnL = dnL / N;


%% function for making the initial parameter vector %%%%%%%%%%%%%%%%%%%%%%%%
function [X] = init_hyp(vbopt, hypinfo)
X = [];
for i=1:length(hypinfo)
  newX = hypinfo(i).hypinvtrans( vbopt.(hypinfo(i).optname) );  
  X = [X; newX(:)];
end

%% function for displaying progress
function stop = outfun(x,optimValues,state)
stop = false;
switch state
  case 'init'
    % do nothing
  case 'iter'
    fprintf('\niter %d: LL=%g; ', optimValues.iteration, -optimValues.fval);
  case 'done'
    fprintf('\nfinal LL=%g; ', -optimValues.fval);
  otherwise
    % do nothing
end



function vbopt = setdefault(vbopt, field, value)
if ~isfield(vbopt, field)
  vbopt.(field) = value;
end

