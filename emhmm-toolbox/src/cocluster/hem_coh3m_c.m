function h3m_out = hem_coh3m_c(h3m_bs,mopt,emopt) 
% H3M Toolbox  (modified function for co-clustering)
% Copyright 2012, Emanuele Coviello, CALab, UCSD
% input: h3m_bs -> list of H3Ms stimuli number * H3M, mopt -> init mode,
%   etc. emopt -> Kb, stimuli number etc.
% output: H3M2 format, in which M*N HMMs
% Notice that the output H3M (H3M2) format is different from original input H3M
% format
% 2016-12-09: ABC - added initialization using splitting
% 2017-05-18: ABC - added support for 'gmm' initialization
% 2017-05-18: ABC - added support for 'gmmNew' initialization
% 2017-08-09: ABC - added 'keep_suboptimal_hmms' option
% 2017-08-17: ABC - added 'full' covariance matrices
% 2020-03-10: ABC - fix warning caused by parfor loop (h3m_new may be uninitialied)
% 2020-03-12: ABC - change 'auto' mode to use hem-ind, gmmNew, and gmmNew2
% 2020-03-12: v0.77 - initial version
% 2020-08-06: v0.80 - fix bug w/ empty HMMs
% 2021-03   : LH - fix for VBHEM

% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2019-07-13
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

% need some functions from the gmm package... 

VERBOSE_MODE = mopt.verbose;
%mopt.K -> group number
if (mopt.K == emopt.subjectsize) || (mopt.K == 0)
    h3m_out = h3m_bs;
    return
end

% ABC 2020-08-06: BUG FIX v0.80: remove empty HMMs for performing initialization
% this removes the association between subjects between stimuli, but is
% fine for initialization.
h3m_bs_no_empty = h3m_bs;
validinds = {};
fprintf('init check: ignore missing hmms: ');
for i=1:emopt.stimulisize
  myinds = find(h3m_bs(i).omega > 0);
  validinds{i} = myinds;
  fprintf(' %d', emopt.subjectsize-length(myinds));
  
  if (length(myinds) < emopt.subjectsize)
    h3m_bs_no_empty(i).K = length(myinds);
    h3m_bs_no_empty(i).omega = h3m_bs_no_empty(i).omega(myinds);
    h3m_bs_no_empty(i).omega = h3m_bs_no_empty(i).omega / sum(h3m_bs_no_empty(i).omega);
    h3m_bs_no_empty(i).hmm = h3m_bs_no_empty(i).hmm(myinds);
  end
end
fprintf('\n');

%{
for i=1:length(h3m_bs_no_empty(3).hmm)
  if length(h3m_bs_no_empty(3).hmm{i}.prior)==1
    h3m_bs_no_empty(3).hmm{i}.emit{1}
  end
end
for k=1:length(h3m_bs_no_empty)
  for i=1:length(h3m_bs_no_empty(k).hmm)
    for j=1:length(h3m_bs_no_empty(k).hmm{i}.prior)
      if all(h3m_bs_no_empty(k).hmm{i}.emit{j}.centres == 0)
        %h3m_bs_no_empty(k).hmm{i}.emit{j}
        fprintf('stim %d, subj %d, emit %d [%d]\n', ...
          k, i, j, h3m_bs_no_empty(k).omega(i)==0)
        
      end
    end
  end
end
%}

%emopt only has one field -> trial=100, inherited from hemopt

LL_best = -inf;

re_do = 0;

% handle different initialization methods
switch (mopt.initmode)
  %% initialize using random base components, or GMM clustering
  case {'base', 'baseem', 'gmm', 'gmmNew', 'gmmNew2', 'kmeans'}
    % typically not executed
    %if (mopt.keep_suboptimal_hmms)
    all_h3ms = cell(1,emopt.trials);
    all_LL   = nan*zeros(1,emopt.trials);
    %end
              
    if (VERBOSE_MODE == 1)
      fprintf('VHEM Trial: ');
    end
    
    parfor t = 1 : emopt.trials
      
      if ~isempty(mopt.seed)
        %fprintf('+ seed %d\n', mopt.seed+t);
        rng(mopt.seed+t, 'twister');
      end
      
      % initialize h3m
      Kb = emopt.subjectsize;
      %            new_idx = randperm(Kb);
      %            new_idx = sort(new_idx);
      
      %            h3m_b1.omega = h3m_b.omega(new_idx);
      %            h3m_b1.hmm = h3m_b.hmm(new_idx);
      %            h3m_b1.K = h3m_b.K;
      
      % ABC 2020-08-06: ignore empty HMMs
      h3m_bs_pointer = h3m_bs_no_empty;
      %mopt.start_time = tic;
      
      h3ms = [];
      
      
      try
        % use less data (400) for initialization
        Kconstant = 400;
        h3ms = [];
        moptpointer = mopt;
        moptpointer.keep_suboptimal_hmms = 0;
        if Kb > Kconstant
          inds = randi(Kb-Kconstant);
          new_idx = inds:(inds+Kconstant-1);
          %initialize all stimuli
          for i=1:emopt.stimulisize
            
            h3m_bs_pointer(i).K = Kconstant;
            h3m_bs_pointer(i).omega = h3m_bs_pointer(i).omega(new_idx);
            h3m_bs_pointer(i).omega = h3m_bs_pointer(i).omega / sum(h3m_bs_pointer(i).omega);
            h3m_bs_pointer(i).hmm = h3m_bs_pointer(i).hmm(new_idx);
            moptpointer.N = mopt.N(i);
            h3m = initialize_hem_h3m_c(h3m_bs_pointer(i),moptpointer);
            h3ms = [h3ms; h3m];
          end
          
        else
          %initialize all stimuli
          for i=1:emopt.stimulisize
              
              if length(mopt.N) > 1
                  moptpointer.N = mopt.N(i);
              else
                  moptpointer.N = mopt.N;
              end
                  
            h3m = initialize_hem_h3m_c(h3m_bs_pointer(i),moptpointer);
            h3ms = [h3ms; h3m];
          end
          
        end
      catch me
        % there is some bug in arkmeans -- need to change it
        % anyways, invalidate this trial
        h3ms = -1;
      end
   
      
      % initialize h3m_new (for parfor loop)
      h3m_new = struct;
      h3m_new.LogL = -inf;            
      
      if isnumeric(h3ms) && (h3ms == -1)        
        % bad initialization, so skip this HEM run
        h3m_new = struct;
        h3m_new.LogL = -inf;
            
      else      
        
        %% run HEM using h3m as the initialization
        % the new reduced h3m is two-dimentional and different from the
        % orginal H3M format
        try
          h3m_new = hem_coh3m_c_step(h3ms,h3m_bs,mopt,emopt);
        catch me
          h3m_new = struct;
          h3m_new.LogL = -inf;
        end
      end
      
      if (VERBOSE_MODE >= 2)
        fprintf('Trial %d\t - loglikelihood: %d\n',t,h3m_new.LogL)
      elseif (VERBOSE_MODE == 1)
        %fprintf('%d ', t);
        %fprintf('.');
        fprintf('trial=%d; LL=%g\n', t, h3m_new.LogL);
        %if mod(t,10)==0
        %  fprintf('%d', t);
        %end
      end
      
      all_h3ms{t} = h3m_new;
      all_LL(t)   = h3m_new.LogL;
    end
    
  case {'hem-ind'}
    %% ABC - still testing
    
    % setup HEM run for each stimuli
    moptind = mopt;
    moptind.termvalue = mopt.termvalue/10;
    moptind.trials    = 1;  %mopt.trials/10;
    moptind.initmode = 'baseem';
    moptind.initopt.mode = 'u';
    moptind.computeStats = 0;
    moptind.verbose = 0;
    moptind.keep_suboptimal_hmms = 0;
    emoptind = emopt;
    emoptind.trials = 1; %emopt.trials/10;
    
    all_h3ms = cell(1,emopt.trials);
    all_LL   = nan*zeros(1,emopt.trials);
    
    % try several times
    parfor t=1:mopt.trials
      
      if ~isempty(mopt.seed)
        rng(mopt.seed+t, 'twister');
        %fprintf('+ set seed %d\n', mopt.seed+t);
      end
      
      mymoptind = moptind;
      
      % run HEM on each stimuli independently first
      indh3m = cell(1,emopt.stimulisize);
      fprintf('init: ');
      for i=1:emopt.stimulisize
        fprintf('.');  %, i, emopt.stimulisize);
        mymoptind.N = mopt.N(i);
        indh3m{i} = hem_h3m_c(h3m_bs(i),mymoptind,emoptind);
      end
      %fprintf('\n');
      
      % match clusters across stimuli
      % Z_combo [stimuli x base x reduced]
      Z_combo = zeros(emopt.stimulisize, emopt.subjectsize, mopt.K);
      for i=1:emopt.stimulisize
        Z_combo(i,:,:) = log(indh3m{i}.Z+1e-30);
      end
      
      newinds = match_clusters(Z_combo, 'k');
      %fprintf('%d ', newinds(:,1)');
      
      % build initialization
      h3ms = [];
      omegall = zeros(1, mopt.K);
      for i=1:emopt.stimulisize
        h3m = struct;
        h3m.K   = indh3m{i}.K;
        h3m.hmm = indh3m{i}.hmm(newinds(i,:));
        h3m.omega = indh3m{i}.omega(newinds(i,:));
        h3m.LogL = -inf;
        h3ms = [h3ms; h3m];
        omegall = omegall + h3m.omega;
      end
      
      % average omega so it is the same across all stimuli
      omegall = omegall / emopt.stimulisize;
      for i=1:emopt.stimulisize
        h3ms(i).omega = omegall;
      end
            
      % run h3m
      h3m_new = hem_coh3m_c_step(h3ms,h3m_bs,mopt,emopt);
      
      if (VERBOSE_MODE >= 2)
        fprintf('HEM-ind\t - loglikelihood: %d\n',h3m_new.LogL)
      elseif (VERBOSE_MODE == 1)
        %fprintf('.');
        fprintf('trial=%d; LL=%g\n', t, h3m_new.LogL);
      end
      
      all_h3ms{t} = h3m_new;
      all_LL(t) = h3m_new.LogL;
    end
    
    
  case {'hem-ind-old'}
    error('not updated to use all_h3ms')
    %% ABC
    % run HEM on each stimuli independently first
    moptind = mopt;
    moptind.termvalue = mopt.termvalue/10;
    moptind.trials    = mopt.trials;
    moptind.initmode = 'baseem';
    moptind.initopt.mode = 'u';
    moptind.computeStats = 0;
    emoptind = emopt;
    emoptind.trials = emopt.trials;
    
    parfor i=1:emopt.stimulisize
      fprintf('=== init %d/%d ===\n', i, emopt.stimulisize);
      if ~isempty(mopt.seed)
        rng(mopt.seed, 'twister');
        %fprintf('+ set seed %d\n', mopt.seed);
      end
      mymoptind = moptind;
      mymoptind.N = mopt.N(i);
      indh3m{i} = hem_h3m_c(h3m_bs(i),mymoptind,emoptind);
    end
    
    % match clusters across stimuli
    % Z_combo [stimuli x base x reduced]
    Z_combo = zeros(emopt.stimulisize, emopt.subjectsize, mopt.K);
    for i=1:emopt.stimulisize
      Z_combo(i,:,:) = log(indh3m{i}.Z+1e-30);
    end
    
    % try different matches (trials)
    all_newinds = {};
    
    %for t=1:emopt.trials
    %newinds = match_clusters(Z_combo, 'r', 0, 0);
    
    %for t=1:emopt.stimulisize
    %  newinds = match_clusters(Z_combo, 'g', t, 0);  % greedy method
    
    %params = logspace(-0.5,0.5,9);
    %params = 1;
    %for t=1:length(params)
    %  newinds = match_clusters(Z_combo, 's', params(t), 1);
    
    for t=1:emopt.trials
      newinds = match_clusters(Z_combo, 'k');
      
      %% this code is same for differnt for loops
      
      % check if this set of inds has been tried
      new_one = 1;
      for q=1:length(all_newinds)
        if all(newinds(:) == all_newinds{q}(:))
          new_one = 0;
          break;
        end
      end
      all_newinds{end+1} = newinds;
      
      if ~new_one
        % check if already tried
        fprintf('skip\n');
        
      else
        fprintf('%d ', newinds(:,1)');
        
        % build initialization
        h3ms = [];
        omegall = zeros(1, mopt.K);
        for i=1:emopt.stimulisize
          h3m.K   = indh3m{i}.K;
          h3m.hmm = indh3m{i}.hmm(newinds(i,:));
          h3m.omega = indh3m{i}.omega(newinds(i,:));
          h3m.LogL = -inf;
          h3ms = [h3ms; h3m];
          omegall = omegall + h3m.omega;
        end
        
        % average omega so it is the same across all stimuli
        omegall = omegall / emopt.stimulisize;
        for i=1:emopt.stimulisize
          h3ms(i).omega = omegall;
        end
        
        % run h3m
        h3m_new = hem_coh3m_c_step(h3ms,h3m_bs,mopt,emopt);
        
        if (VERBOSE_MODE >= 2)
          fprintf('HEM-ind\t - loglikelihood: %d\n',h3m_new.LogL)
        elseif (VERBOSE_MODE == 1)
          %fprintf('.');
          fprintf('LL=%g\n', h3m_new.LogL);
        end
        
        % save if the new model has the best LL
        if (t == 1) || (h3m_new.LogL > LL_best)
          LL_best = h3m_new.LogL;
          h3m_out = h3m_new;
          t_best = t;
        end
        
      end
      
    end
    fprintf('Best run is trial %d: loglikelihood: %d\n', t_best, LL_best)
    
  otherwise
    error('unknown initialization');
end

%% keep subooptimal HMMs
if (mopt.keep_suboptimal_hmms)
  diffthresh = 2*mopt.termvalue*10;
  unique_inds = uniqueLL(all_LL, diffthresh, 1);
  unique_h3ms = {all_h3ms{unique_inds}};
  suboptimal_h3ms = unique_h3ms;
  trials_LL = all_LL;
  suboptimal_h3ms_inds = unique_inds;
end


%% check for degenerate solution
%% if trials is 1,there will be a bug
if mopt.trials>1
[all_LL, all_h3ms] = check_degenerate(all_LL, all_h3ms);
end
    
%% get the best model
[~, t_best] = max(all_LL);
LL_best = all_LL(t_best);
h3m_out = all_h3ms{t_best};

%% lh: keep the top h3ms as initialzation for vbhem-cocluster
if (isfield(mopt, 'keepTopK') && mopt.keepTopK)
    TOP = min(mopt.keepTopK, mopt.trials);
    [~, inds] = sort(all_LL, 'descend');
    topK_h3ms = all_h3ms(inds(1:TOP));
    topK_h3ms_inds = inds(1:TOP);
    trials_LL2 = all_LL;
    h3m_out.topK_h3ms = topK_h3ms;
    h3m_out.suboptimal_h3ms_inds = topK_h3ms_inds;
    h3m_out.trials_LL = trials_LL2;
end

if (VERBOSE_MODE >= 1)
  fprintf('\nBest run is %d: LL=%g\n',t_best,LL_best)
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [g1, g2] = split_gauss(g, f)
% return 2 new gaussians, split from g
if nargin<2
  f = 1;
end

g1 = g;
g2 = g;

cv = g.covars;
switch(g.covar_type)
  case 'diag'
    [~,mi] = max(cv);
    delta_mn = zeros(size(cv));
    delta_mn(mi) = sqrt(cv(mi));
    new_cv = cv;
    new_cv(mi) = new_cv(mi) / ((2*f)^2);
  case 'full'
    error('not supported yet');
end

g1.centres = g1.centres + f*delta_mn;
g2.centres = g2.centres - f*delta_mn;
g1.covars = new_cv;
g2.covars = new_cv;


function A = soften_A(A, f)

A = A + f;
A = bsxfun(@times, A, 1./sum(A,2));


%% find outliers (degenerate solutions)
function [all_LL, all_h3ms] = check_degenerate(all_LL, all_h3ms)

outliers = find_degenerate_outliers(all_LL, all_h3ms);

if length(outliers)>0
  figure
  tmp = randn(size(all_LL));
  plot(all_LL, tmp, 'bx');
  hold on
  plot(all_LL(outliers), tmp(outliers), 'rx');
  hold off
  drawnow
  
  fprintf('*** removing outliers: ***\n');
  outliers
%  keyboard
  all_LL(outliers) = -inf;
end

