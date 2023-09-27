function [group_hmms] = vhem_cluster(hmms, K, S, hemopt)
% vhem_cluster - cluster hmms into groups using VHEM algorithm
%
% [group_hmms] = vhem_cluster(hmms, K, S, hemopt)
%
% This function uses a modified version of the H3M toolbox.
%
% INPUTS
%      hmms = Nx1 cell array of HMMs, learned with vbhmm_learn.
%
%         K = number of groups to use (required)
%         S = number of states (ROIs) in a group HMM 
%             (default=[], use the median number of states in input HMMs)
%
%    hemopt = structure containing other options as below:
%      hemopt.trials  = number of trials with random initialization to run (default=100)
%      hemopt.tau     = the length of the data sequences (default=10)
%      hemopt.sortclusters = '' (no sorting)
%                          = 'p' - sort ROIs by prior frequency
%                          = 'f' - sort ROIs by most likely fixation path [default]
%                            (see vbhmm_standardize for more options)
%      hemopt.remove_empty = 1 = remove empty ROIs in the input HMMs before running VHEM [default]
%                            0 = use the original input HMMs
%      hemopt.verbose  = 0 - no messages
%                        1 - progress messages [default]
%                        2 - more messages
%                        3 - debugging
%      hemopt.seed    = seed state for random number generator (required for reproducible results)
%
%    More VHEM options:
%      hemopt.initmode = 'auto' - try multiple initialization methods and select the best [default]
%                        'baseem' - randomly select base HMM emissions as initialization
%                        'gmmNew' - estimate GMM from the HMM emissions, and initialize each group with the same GMM.
%                        'gmmNew2' - esitmate a larger GMM from the HMM emissions, and initialize each group with different GMMs
%                        'base'  - randomly select base HMMs components as initialization 
%      hemopt.initopt.iter   = number of iterations for GMM initialization (default=30)
%      hemopt.initopt.trials = number of trials for GMM initialization (default=4)
%      hemopt.reg_cov = regularization for covariance matrix (default=0.001)
%      hemopt.max_iter = maximum number of iterations
%
% OUTPUT
%   group_hmm = structure containing the group HMMs and information
%   group_hmm.hmms       = cell array of group HMMs {1xK}
%   group_hmm.label      = the cluster assignment for each input HMM [1xN]
%   group_hmm.groups     = indices of input HMMs in each group {1xK}
%   group_hmm.group_size = size of each group [1xK]
%   group_hmm.Z          = probabilities of each input HMM belonging to a group [Nx2]
%   group_hmm.LogL       = log-likelihood score
%
%   statistics
%   group_hmm.hmms{j}.stats.
%      emit_vcounts    = number of virtual samples assigned to each group hmm ROI (it sums to num virtual samples).
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-01-13
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong

% INTERNAL options
%    hemopt.useMEX       = 1 - use MEX implementations [default]
%                        = 0 - do not use MEX
%   hemopt.keep_suboptimal_hmms = 0 - do not keep
%                               = 1 - keep in vb_hmm.suboptimal_hmms{i}
%   hemopt.computeStats = 1 - compute stats for MANOVA 
%                       = 0 - don't compute [default]
%    this option will generate STATS FOR MANOVA:
%      emit{i}.weights = the percentage of input HMM ROIs that are assigned to each group hmm ROI
%      emit{i}.mu2     = second moment of the input HMM ROI means, for each group hmm ROI.
%      emit{i}.mu      = first moment of the input HMM ROI means
%      emit{i}.NROIs   = the number of individual ROIs contributing to each group hmm ROI


% 2016-11-21: ABC - initial version (modified from Tim)
% 2016-12-09: ABC - added initialization using splitting
% 2016-01-10: ABC - added initialization using base emissions (baseem)
% 2017-01-11: ABC - added verbose option
% 2017-01-19: ABC - added fixation duration
% 2017-05-18: ABC - added support for 'gmm' initialization -- perhaps not very good,
%                   since each HMM component has the same Gausian.
%                 - added support for 'gmmNew' initialization - pool the emissions, and
%                   estimate a GMM. Each comp of the GMM is used as an emission.
% 2017-08-08: ABC - use MEX implementation to speed up E-step
% 2017-08-09: ABC - added "keep_suboptimal_hmms" option
% 2017-08-12: ABC - added "auto" option for initmode
% 2017-08-14: ABC - added support for full covariance matrix [now default]
% 2017-11-25: ABC - small bug fix when checking covar_type
% 2018-02-13: v0.72 - handle degenerate HMM states
% 2018-05-24: v0.73 - minor bug fix on Zomega
% 2018-05-25: v0.73 - add stats for MANOVA test.
% 2018-06-08: v0.73 - bug fix in VHEM when one ROI.
% 2018-06-13: v0.73 - bug fix in VHEM when tau=1
% 2018-06-13: v0.73 - option to remove empty ROIs in input HMMs [default to on]
% 2019-01-12: v0.75 - added hemopt.seed option
% 2019-02-10: v0.75 - use parfor to improve speed
% 2019-02-20: v0.75 - add version number
% 2021-07-22: v0.78 - support for fixed ROIs
% 2021-09-10: v0.80 - add matlab version number in output struct

if nargin<3
  S = [];
end
if nargin<4
  hemopt = struct;
end

% number of clusters
hemopt.K = K;
hemopt = setdefault(hemopt, 'verbose', 1);
VERBOSE_MODE = hemopt.verbose;


% v0.78 - handle fixed ROIs
if ~isempty(hmms{1}.vbopt.fixed_rois)
  if VERBOSE_MODE >= 1
    fprintf('input HMMs have fixed ROIs, using fixed ROI seting\n');
  end
  if isfield(hemopt, 'remove_empty') && (hemopt.remove_empty == 1)
    error('fixed ROIs: remove_empty is set to 1, but should not be set');
  end
  hemopt.remove_empty = 0;
  
  if isfield(hemopt, 'sortclusters') && ~isempty(hemopt.sortclusters)
    error('fixed ROIs: sortclusters is set, but it is not suggested to set since it will reorder the ROIs');
  end
  hemopt.sortclusters = '';
  
  if isfield(hemopt, 'initmode') && ~isempty(hemopt.initmode)
    error('fixed ROIs: initmode is set, but it is not be set');
  end
  hemopt.initmode = 'base';
end


hemopt = setdefault(hemopt, 'remove_empty', 1);  % remove empty ROIs


% v0.73 - remove empty states
% do this before setting number of states
% v0.74 - tidy up
if (hemopt.remove_empty)
  if (hemopt.verbose)
    fprintf('Checking input HMMS: ');
  end
  for i=1:length(hmms)
    [hmms{i}, zi] = vbhmm_remove_empty(hmms{i}, 0, 1e-3);
    if ~isempty(zi)
      if (hemopt.verbose)
        fprintf('%d: removed states', i);
        fprintf(' %d', zi);
        fprintf('; ');
      end
    end
    
  end
  if (hemopt.verbose)
    fprintf('done\n');
  end
end

% number of states in a group HMM
if isempty(S)
  % get median number
  hmmsS = cellfun(@(x) length(x.prior), hmms);
  S = ceil(median(hmmsS)); % if in the middle, then take the larger S
  
  if (VERBOSE_MODE >= 1)
    fprintf('using median number of states: %d\n', S);
  end
end
hemopt.N = S;
%hemopt = setdefault(hemopt, 'N',  3); % Number of states of each HMM

% setable parameters (could be changed)
hemopt = setdefault(hemopt, 'trials', 100);     % number of trials to run
hemopt = setdefault(hemopt, 'reg_cov', 0.001);   % covariance regularizer
hemopt = setdefault(hemopt, 'termmode', 'L');   % rule to terminate EM iterations
hemopt = setdefault(hemopt, 'termvalue', 1e-5);
hemopt = setdefault(hemopt, 'max_iter', 100);    % max number of iterations
hemopt = setdefault(hemopt, 'min_iter', 1);     % min number of iterations
hemopt = setdefault(hemopt, 'sortclusters', 'f');
hemopt = setdefault(hemopt, 'initmode', 'auto');  % initialization mode
hemopt = setdefault(hemopt, 'tau', 10); % temporal length of virtual samples
hemopt = setdefault(hemopt, 'seed', []); % rng seed


% standard parameters (not usually changed)
hemopt = setdefault(hemopt, 'Nv', 100); % number of virtual samples
hemopt = setdefault(hemopt, 'initopt', struct);   % options for 'gmm' initaliation
hemopt.initopt = setdefault(hemopt.initopt, 'iter', 30);  % number of GMM iterations for init
hemopt.initopt = setdefault(hemopt.initopt, 'trials', 4); % number of GMM trials for init
hemopt.initopt = setdefault(hemopt.initopt, 'mode', ''); % other options
if isempty(hemopt.initopt.mode)
  % set defaults for each initialization type
  switch(hemopt.initmode)
    case 'baseem'
      hemopt.initopt.mode = 'u';
    case 'gmmNew'
      hemopt.initopt.mode = 'r0';
    case 'gmmNew2'
      hemopt.initopt.mode = 'u0';  
  end
end
hemopt = setdefault(hemopt, 'inf_norm', 'nt');   % normalization before calculating Z ('nt'=tau*N/K). This makes the probabilites less peaky (1 or 0).
hemopt = setdefault(hemopt, 'smooth', 1);        % smoothing parameter - for expected log-likelihood
hemopt = setdefault(hemopt, 'useMEX', 1);  % use MEX implementation
hemopt = setdefault(hemopt, 'keep_suboptimal_hmms', 0); % keep other suboptimal solutions
hemopt = setdefault(hemopt, 'computeStats', 0); % compute stats for MANOVA

% fixed parameters (not usually changed)
hemopt.emit.type = 'gmm';
hemopt.emit = setdefault(hemopt.emit, 'covar_type', 'full');
hemopt.M = 1; % number of mixture components in each GMM for an emission (should always be 1)

% a separate structure
emopt.trials = hemopt.trials;

% v0.75 - set seed
if isempty(hemopt.seed)
  error('hemopt.seed needs to be set to make reproducible results.');
end
if ~isempty(hemopt.seed)
  if (hemopt.verbose)
    fprintf('+ set seed to %d\n', hemopt.seed);
  end
  rng(hemopt.seed, 'twister');
end

% convert list of HMMs into an H3M
H3M = hmms_to_h3m(hmms, hemopt.emit.covar_type);

if ~strcmp(hemopt.initmode, 'auto')
  %% run VHEM clustering
  h3m_out = hem_h3m_c(H3M, hemopt, emopt);
  
else
  %% auto initialization - try these methods
  initmodes = {'baseem', 'gmmNew', 'gmmNew2'};
  initopts = {'u', 'r0', 'u0'};
  
  for ii=1:length(initmodes)
    hemopt2 = hemopt;
    hemopt2.initmode = initmodes{ii};
    hemopt2.initopt.mode = initopts{ii};
    if (hemopt.verbose)
      fprintf('auto initialization: trying %s: ', hemopt2.initmode);
    end
    h3m_out_trials{ii} = hem_h3m_c(H3M, hemopt2, emopt);
    trials_LL(ii) = h3m_out_trials{ii}.LogL;
  end
  
  % get best method
  [maxLL, ind] = max(trials_LL);
  
  if (hemopt.verbose)
    fprintf('best init was %s; LL=%g\n', initmodes{ind}, trials_LL(ind));
  end
  h3m_out = h3m_out_trials{ind};
end

% convert back to our format
group_hmms = h3m_to_hmms(h3m_out);


%% v0.73 - collect more stats for MANOVA test
if (hemopt.computeStats)
  warning('TODO: make the stats permute in vbhmm_standardize');
  
  % get total number of ROIs
  totIndROIs = 0;
  for i=1:length(hmms)
    totIndROIs = totIndROIs + size(hmms{i}.trans, 1);
  end
  group_hmms.stats.totIndROIs = totIndROIs;
  
  % collect all counts
  counts = [];
  for j=1:K
    counts = [counts group_hmms.hmms{j}.stats.emit_vcounts];
  end
  Tcounts = sum(counts);
  
  for j=1:K
    for i=1:length(group_hmms.hmms{j}.stats.emit_vcounts)
      % normalize to get weights
      group_hmms.hmms{j}.stats.emit{i}.weights = group_hmms.hmms{j}.stats.emit_vcounts(i) / Tcounts;
      
      % get Number of ROIs
      group_hmms.hmms{j}.stats.emit{i}.NROIs = totIndROIs*group_hmms.hmms{j}.stats.emit{i}.weights;
    end
  end
end

% sort clusters
if ~isempty(hemopt.sortclusters)
  group_hmms = vbhmm_standardize(group_hmms, hemopt.sortclusters);
end
  
% keep suboptimal solutions
if (hemopt.keep_suboptimal_hmms)
  group_hmms.suboptimal_hmms = cell(1,length(h3m_out.suboptimal_h3ms));
  for i=1:length(h3m_out.suboptimal_h3ms)
    tmp = h3m_to_hmms(h3m_out.suboptimal_h3ms{i});
    if ~isempty(hemopt.sortclusters)
      tmp = vbhmm_standardize(tmp, hemopt.sortclusters);
    end
    group_hmms.suboptimal_hmms{i} = tmp;
  end
  group_hmms.trials_LL = h3m_out.trials_LL;
end


% save parameters
group_hmms.hemopt = hemopt;

% add version number
group_hmms.emhmm_version = emhmm_version();
group_hmms.matlab_version = version();

function hemopt = setdefault(hemopt, field, value)
if ~isfield(hemopt, field)
  hemopt.(field) = value;
end


