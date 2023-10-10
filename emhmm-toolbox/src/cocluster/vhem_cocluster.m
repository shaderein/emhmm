function [group_hmms, h3m_out_trials] = vhem_cocluster(hmms, K, S, hemopt)
% vhem_cocluster - co-clustering hmms into groups using VHEM algorithm
%
% [group_hmms] = vhem_cocluster(hmms, K, S, hemopt)
%
% Perform co-clustering on sets of HMMs
%
% INPUTS
%      hmms = [Stimuli(Trial) x Subject] cell array of HMMs, learned with vbhmm_learn.
%         K = number of groups to use (required)
%         S = number of states (ROIs) in a group HMM 
%             (default=[], use the median number of states in input HMMs)
%
%    hemopt = structure containing other options as below:
%      hemopt.trials  = number of trials with random initialization to run (default=100)
%      hemopt.tau     = the length of the data sequences (default=10)
%      hemopt.sortclusters = '' (no sorting, default)
%                          = 'p' - sort ROIs by prior frequency
%                          = 'f' - sort ROIs by most likely fixation path [default]
%                            (see vbhmm_standardize for more options)
%      hemopt.verbose  = 0 - no messages
%                        1 - progress messages [default]
%                        2 - more messages
%                        3 - debugging
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
%   group_hmm{i} = cell array containing the group HMMs and information for i-th stimuli
%   group_hmm{i}.hmms       = cell array of HMMs (K)
%   group_hmm{i}.label      = the cluster assignment for each input HMM [1xN]
%   group_hmm{i}.groups     = indices of input HMMs in each group (K)
%   group_hmm{i}.group_size = size of each group (K)
%   group_hmm{i}.Z          = probabilities of each input HMM belonging to a group [Nx2]
%   group_hmm{i}.LogL       = log-likelihood score (1)
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2019-07-13
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

% INTERNAL options
%    hemopt.useMEX       = 1 - use MEX implementations [default]
%                        = 0 - do not use MEX
%   hemopt.keep_suboptimal_hmms = 0 - do not keep
%                               = 1 - keep in vb_hmm.suboptimal_hmms{i}

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
% 2018-11-23: ABC v0.74 - handle empty HMMs (hmm2h3m, initialization).
% 2020-03-12:     v0.77 - initial version for co-clustering
% 2020-08-06:     v0.80 - fix bug causing empty ROIs when many trials are missing
% 2020-09-07:     v0.80 - fix bug for initialization with some empty HMMs (missing trials)
% 2021-09-10:           - add version numbers

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

% ABC - added %%%%%%%%%%%%
hemopt = setdefault(hemopt, 'remove_empty', 1);  % remove empty ROIs
% v0.73 - remove empty states
% do this before setting number of states
if (hemopt.remove_empty)
  fprintf('Checking input HMMS: ');
  for i=1:size(hmms,1)
    for j=1:size(hmms,2)
      if ~isempty(hmms{i,j})
        [hmms{i,j}, zi] = vbhmm_remove_empty(hmms{i,j}, 0, 1e-3);
        if ~isempty(zi)
          fprintf('(%d,%d) removed states', i,j);
          fprintf(' %d', zi);
          fprintf('; \n');
        end
      end
    end
  end
  fprintf(' done\n');
end
%%%%%%%%%%%%%%%%

% number of states in a group HMM
sizeall = size(hmms);
stimulisize = sizeall(1);
subjectsize = sizeall(2);

% for each stimuli, check for empty HMMs
hmms_validmask = zeros(stimulisize, subjectsize);
hmms_validinds = {};
for i=1:stimulisize
  tmp = ~cellfun(@isempty, hmms(i,:));
  hmms_validmask(i,:) = tmp;
  hmms_validinds{i} = find(tmp);
end

hemopt.N = [];
if isempty(S)
    % iterate through all stimuli (trials)
    for i=1:stimulisize 
      v = hmms_validinds{i};
      % get median number
      hmmsS = cellfun(@(x) length(x.prior), hmms(i,v));
      S = ceil(median(hmmsS)); % if in the middle, then take the larger S

      if (VERBOSE_MODE > 1)
        fprintf('using median number of states: %d\n', S);
      end

      hemopt.N = [hemopt.N; S];
    end
else
    hemopt.N = S;
end
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

% fixed parameters (not usually changed)
hemopt.emit.type = 'gmm';
hemopt.emit = setdefault(hemopt.emit, 'covar_type', 'full');
hemopt.M = 1; % number of mixture components in each GMM for an emission (should always be 1)


% a separate structure
emopt.trials = hemopt.trials;

% convert list of HMMs into an H3M
H3Ms = [];
for i=1:stimulisize
    H3M = hmms_to_h3m(hmms(i,:), hemopt.emit.covar_type);
    H3Ms = [H3Ms; H3M];
end

emopt.stimulisize = stimulisize;
emopt.subjectsize = subjectsize;

if ~strcmp(hemopt.initmode, 'auto')
  %% run VHEM clustering
  h3m_out = hem_coh3m_c(H3Ms, hemopt, emopt);
  h3m_out_trials = {h3m_out};
  
else
  %% auto initialization - try these methods
  % original VHEM
  %initmodes = {'baseem', 'gmmNew', 'gmmNew2'};
  %initopts = {'u','r0', 'u0'};
  
  % testing hem-ind
  %initmodes = {'baseem','gmmNew2'};
  %initopts = {'u','u0'};
  %initmodes = {'hem-ind', 'baseem'};
  %initopts = {'', 'u'};
  
  % try all methods (still testing  -- baseem is not stable)
  %initmodes = {'hem-ind', 'baseem', 'gmmNew', 'gmmNew2'};
  %initopts = {'', 'u','r0', 'u0'};
  
  initmodes = {'hem-ind', 'gmmNew', 'gmmNew2'};
  initopts = {'','r0', 'u0'};
  
  
  %warning('debugging');
  %initmodes = {'gmmNew', 'gmmNew2'};
  %initopts = {'r0', 'u0'};
 
  for ii=1:length(initmodes)
    hemopt2 = hemopt;
    hemopt2.initmode = initmodes{ii};
    hemopt2.initopt.mode = initopts{ii};
    fprintf('auto initialization: trying %s: ', hemopt2.initmode);
    h3m_out_trials{ii} = hem_coh3m_c(H3Ms, hemopt2, emopt);
    trials_LL(ii) = h3m_out_trials{ii}.LogL;
  end
  
  % get best method
  [maxLL, ind] = max(trials_LL);
  
  fprintf('best init was %s; LL=%g\n', initmodes{ind}, trials_LL(ind));
  h3m_out = h3m_out_trials{ind};
end

group_hmms = cell(stimulisize,1);
for i=1:stimulisize
    %make it compatible with LogL and LogLs
    h3m_out.h3m_rs(i).LogLs = h3m_out.h3m_rs(i).LogL;
    grouphmm = h3m_to_hmms(h3m_out.h3m_rs(i));
    
    % add version number
    grouphmm.emhmm_version = emhmm_version();
    grouphmm.matlab_version = version();
    
    group_hmms{i} = grouphmm;
end

% sort clusters
if ~isempty(hemopt.sortclusters)
  for i=1:stimulisize
    group_hmms{i} = vbhmm_standardize(group_hmms{i}, hemopt.sortclusters);
    group_hmms{i}.hemopt = hemopt;
  end
end

fprintf('Groups:\n');

for i=1:K
  fprintf('  Group %d: %s\n', i, num2str(group_hmms{1}.groups{i}));
end

function hemopt = setdefault(hemopt, field, value)
if ~isfield(hemopt, field)
  hemopt.(field) = value;
end