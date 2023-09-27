function [h3m_r]= vbhem_cluster(hmms, K, S, vbhemopt)
% vbhem_cluster - cluster HMMs using variational Bayesian HEM
%
%  [group_vbhmms]= vbhem_cluster(hmms, K, S, vbhemopt)
%  
% INPUTS
%      hmms = Nx1 cell array of HMMs, learned with vbhmm_learn.
%
%         K = number or vecter of number clusters  (required)  
%         S = number or vecter of states (ROIs) in one HMM     
%       if a vector: automatically selects the number of hidden states from the given values. 
%                    The model (K,S) with highest log-likelihood is selected.
%             
%  vbhembopt = structure containing other options as below: (similar with vbhmm_learn)
%
%   VB hyper-parameters: each hmm has the same initial hyper-parameters
%     vbhemopt.alpha0 = Dirichlet distribution concentration parameter on cluster probabilities (default=1)
%                   -- large value encourages uniform prior, small value encourages concentrated prior 
%                   Another way to think of it is in terms of "virtual" samples -- A typical way to 
%                   estimate the probability of something is to count the number of samples that it
%                   occurs and then divide by the total number of samples, 
%                   i.e. # times it occurred / # samples.  
%                   The alpha parameter of the Dirichlet adds a virtual sample to this estimate, 
%                   so that the probability estimate is (# times it occurred + alpha) / # samples.  
%                   So for small alpha, then the model will basically just do what the data says.  
%                   For large alpha, it forces all the probabilities to be very similar (i.e., uniform).
%     vbhemopt.eta0  = Dirichlet distribution concentration parameter on the initial state probabilities (default=1)
%     vbhemopt.epsilon0 = Dirichlet distribution for rows of the transition matrix (default=1). 
%                     The meaning is similar to alpha, but for the transition probabilities.
%     vbhemopt.m0   = prior mean - it should be dimension D;
%                   for D=2 (fixation location), default=[256;192] -- typically it should be at the center of the image.
%                   for D=3 (location & duration), default=[256,192,250]
%     vbhemopt.W0    = size of the inverse Wishart distribution (default=0.005).
%                   If a scalar, W is assumed to be isotropic: W = vbhemopt.W*I.
%                   If a vector, then W is a diagonal matrix:  W = diag(vbhemopt.W);
%                   This determines the covariance of the ROI prior.
%                   For example W=0.005 --> covariance of ROI is 200 --> ROI has standard devation of 14.
%     vbhemopt.v0     = dof of inverse Wishart, v > D-1. (default=5) -- larger values give preference to diagonal covariance matrices.
%     vbhemopt.lambda0  = Wishart concentration parameter -- large value encourages use of
%                   prior mean & W, while small value ignores them. (default=1)
%     vbhemopt.learn_hyps = 1 - estimate these hyperparameters from the data, use the provided values as the initialization (default).
%                           0 - use fixed hyperparameters from vbhemopt.
%
%   VBHEM Algorithm parameters
%     vbhemopt.initmode     = initialization method (default='baseem')
%                             follow the initialization ways in VHEM, using
%                             the parameters in the h3m to update the
%                             hyperparameters.
%                            'random' - randomly set assignment varibles z and phi 
%                            'baseem' - randomly select a set of base emissions
%                            'inith3m' - specify a H3M as the initialization
%                            'gmmNew', 'gmmNew2' - initialize emissions using HEM-G3M estimated
%                            'gmmNew_r' - combine gmmNew and random,  using the information of mapping in the HEM-G3M
%                            'wtkmeans' - use weighted k-means
%     vbhemopt.trials    = number of trials for 'random' initialization (default=100)
%                              
%     vbhemopt.maxIter      = max number of iterations (default=100)
%     vbhemopt.minDiff      = tolerence for convergence (default=1e-5)
%     vbhemopt.showplot     = show plots (default=0)
%     vbhemopt.sortclusters = '' - no sorting 
%                          'e' - sort ROIs by emission probabilites
%                          'p' - sort ROIs by prior probability
%                          'f' - sort ROIs by most-likely fixation path [default]
%                              (see vbhmm_standardize for more options)
%     vbhemopt.prune_hmms   = 1 - prune empty clusters and empty states from the final HMM (default)
%                             0 - don't prune
%     vbhemopt.remove_empty = 1 - remove empty input HMMs before clustering  (default)
%                             0 - don't remove
%     vbhemopt.verbose   = 0 - no messages
%                        = 1 - a few messages showing progress [default]
%                        = 2 - more messages
%                        = 3 - debugging
%
% OUTPUTS
%   group_vbhmms = structure containing the hmms in the reduced model and related information
%   group_vbhmms.omega  = omega [Kr x 1]
%   group_vbhmms.hmms   = the hmms, cell [Kr x 1]
%   group_vbhmms.hmms{j}.prior       = prior probababilies [Sr x 1]
%   group_vbhmms.hmms{j}.trans       = transition matrix [Sr x Sr]: trans(i,j) is the P(x_t=j | x_{t-1} = i)
%   group_vbhmms.hmms{j}.pdf{j}.mean = ROI mean [1 x D]
%   group_vbhmms.hmms{j}.pdf{j}.cov  = covariances [D x D]
%   group_vbhmms.hmms{j}.LL          = log-likelihood of data
%   group_vbhmms.label       = the cluster for each input HMM
%   group_vbhmms.groups      = the list of clusters (cell array 1xK)
%
% the following are related to the search over all K,S pairs
%   group_vbhmms.model_LL    = [NK x NS] - LL for each (K,S) pair
%   group_vbhmms.model_KS    = [NK x NS] - the (K,S) pair for each entry
%   group_vbhmms.model_bestK = the best K selected
%   group_vbhmms.model_bestS = the best S selected
%   group_vbhmms.model_all   = the learned model for each (K,S) pair
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2021-07-29
% Antoni B. Chan, Janet H. Hsiao, Lan Hui
% City University of Hong Kong, University of Hong Kong

% 2021-07-29: v0.80 - Hui - initial version
%                   - ABC - added "remove_empty" option
%                         - search K,S in one big loop.
%                         - add version numbers


% do not know now!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%   h3m_r.hmm{j}.gamma       = responsibilities gamma{i}=[KxT] -- probability of each fixation belonging to an ROI.
%   h3m_r.hmm{j}.M           = transition counts [K x K]
%   h3m_r.hmm{j}.N1          = prior counts [K x 1]
%   h3m_r.hmm{j}.N           = cluster sizes [K x 1] (emission counts)
% OUTPUT
%   group_hmm = structure containing the group HMMs and information
%   group_hmm.hmms       = cell array of group HMMs {1xK}
%   group_hmm.label      = the cluster assignment for each input HMM [1xN]
%   group_hmm.groups     = indices of input HMMs in each group {1xK}
%   group_hmm.group_size = size of each group [1xK]
%   group_hmm.Z          = probabilities of each input HMM belonging to a group [Nx2]
%   group_hmm.LogL       = log-likelihood score
%

if nargin<4
  vbhemopt = struct;
end

%% set default values

% number of clusters
vbhemopt.K = K;

% number of states
vbhemopt.S = S;

vbhemopt = setdefault(vbhemopt, 'verbose', 1);
vbhemopt = setdefault(vbhemopt, 'remove_empty', 1);

if (vbhemopt.remove_empty)
  %if (vbhemopt.verbose)
  %  fprintf('Checking input HMMS: ');
  %end
  for i=1:length(hmms)
    [hmms{i}, zi] = vbhmm_remove_empty(hmms{i}, 0, 1e-3);
    if ~isempty(zi)
      if (vbhemopt.verbose)
        fprintf('input HMM %d: removed states', i);
        fprintf(' %d', zi);
        fprintf('; ');
      end
    end
    
  end
  %if (vbhemopt.verbose)
  %  fprintf('done\n');
  %end
end



% assume the data has same dimensions
D = size(hmms{1}.pdf{1}.mean,2);
switch(D)
  case 2
    defmu = [256;192]; % for face data
    %defmu = zeros(D,1);
  case 3
    defmu = [256;192;150];
  otherwise
    defmu = zeros(D,1);
    %warning(sprintf('no default mu0 for D=%d', D));
end

vbhemopt = setdefault(vbhemopt, 'alpha0', 1);
vbhemopt = setdefault(vbhemopt, 'eta0', 1);
vbhemopt = setdefault(vbhemopt, 'epsilon0', 1);
vbhemopt = setdefault(vbhemopt, 'm0',   defmu); % hyper-parameter for the mean
vbhemopt = setdefault(vbhemopt, 'W0',     0.005); % the inverse of the variance of the dimensions
vbhemopt = setdefault(vbhemopt, 'lambda0',  1);
vbhemopt = setdefault(vbhemopt, 'v0',     5);

% setable parameters (could be changed)
vbhemopt = setdefault(vbhemopt, 'trials', 100);     % number of trials to run  in vhem is 'numtrials'
vbhemopt = setdefault(vbhemopt, 'termmode', 'L');   % rule to terminate EM iterations
vbhemopt = setdefault(vbhemopt, 'termvalue', 1e-5);
vbhemopt = setdefault(vbhemopt, 'max_iter', 200);    % max number of iterations
vbhemopt = setdefault(vbhemopt, 'min_iter', 1);     % min number of iterations
vbhemopt = setdefault(vbhemopt, 'sortclusters', 'f');
vbhemopt = setdefault(vbhemopt, 'initmode',  'auto');  % initialization mode
vbhemopt = setdefault(vbhemopt, 'minDiff',   1e-5);
vbhemopt = setdefault(vbhemopt, 'showplot',  0);
vbhemopt = setdefault(vbhemopt, 'random_gmm_opt', {});


% standard parameters (not usually changed)
vbhemopt = setdefault(vbhemopt, 'Nv', 10); % number of virtual samples
vbhemopt = setdefault(vbhemopt, 'tau', 10); % temporal length of virtual samples
vbhemopt = setdefault(vbhemopt, 'initopt', struct);   % options for 'gmm' initaliation (unused)
vbhemopt.initopt = setdefault(vbhemopt.initopt, 'iter', 30);  % number of GMM iterations for init (unused)
vbhemopt.initopt = setdefault(vbhemopt.initopt, 'trials', 4); % number of GMM trials for init (unused)
vbhemopt.initopt = setdefault(vbhemopt.initopt, 'mode', ''); % other options

vbhemopt = setdefault(vbhemopt, 'learn_hyps', 1);

% minimum/maximum values of hyperparameters, for numerical stability
vbhemopt = setdefault(vbhemopt, 'hyps_max', struct);
hyps_max = vbhemopt.hyps_max;
hyps_max = setdefault(hyps_max, 'alpha0',   1.0686e+13);  % exp(30)
hyps_max = setdefault(hyps_max, 'eta0',   1.0686e+13);  % exp(30)
hyps_max = setdefault(hyps_max, 'epsilon0', 1.0686e+13);  % exp(30)
hyps_max = setdefault(hyps_max, 'v0',       1e4);
hyps_max = setdefault(hyps_max, 'lambda0',    1.0686e+13);  % exp(30)
hyps_max = setdefault(hyps_max, 'W0',       1.0686e+13);  % exp(30)
vbhemopt.hyps_max = hyps_max;

vbhemopt = setdefault(vbhemopt, 'hyps_min', struct);
hyps_min = vbhemopt.hyps_min;
hyps_min = setdefault(hyps_min, 'alpha0',    1.0686e-13);  % exp(-30)
hyps_min = setdefault(hyps_min, 'eta0',    1.0686e-13);  % exp(-30)
hyps_min = setdefault(hyps_min, 'epsilon0', 1.0686e-13);  % exp(-30)
hyps_min = setdefault(hyps_min, 'v0',       2.0612e-09+D-1);  % exp(-20)+dim-1 - should be okay to ~10000 dimensions
hyps_min = setdefault(hyps_min, 'lambda0',    1.0686e-13);  % exp(-30)
hyps_min = setdefault(hyps_min, 'W0',       1.0686e-13);  % exp(-30)
vbhemopt.hyps_min = hyps_min;


vbhemopt = setdefault(vbhemopt, 'calc_LLderiv', 0);
vbhemopt = setdefault(vbhemopt, 'keep_suboptimal_hmms', 0);
vbhemopt = setdefault(vbhemopt, 'minimizer', 'minimize-lbfgs');
vbhemopt = setdefault(vbhemopt, 'canuseMEX', 1);
vbhemopt = setdefault(vbhemopt, 'verbose_prefix', '');
vbhemopt = setdefault(vbhemopt, 'keep_best_random_trial', 1);
vbhemopt = setdefault(vbhemopt, 'seed', []);
vbhemopt = setdefault(vbhemopt, 'use_post', 1);

vbhemopt = setdefault(vbhemopt, 'prune_hmms', 1);

if isempty(vbhemopt.initopt.mode)
  % set defaults for each initialization type
  switch(vbhemopt.initmode)
    case 'baseem'
      vbhemopt.initopt.mode = 'u';
    case {'gmmNew','wtkmeans','con-wtkmeans'}
      vbhemopt.initopt.mode = 'r0';
    case {'gmmNew2','gmmNew2_r'}
      vbhemopt.initopt.mode = 'u0';
  end
end


% fixed parameters (not changed)
vbhemopt.emit.type= 'gmm';
vbhemopt.emit = setdefault(vbhemopt.emit, 'covar_type', 'full');

%vbhemopt.emit.covar_type= 'diag';  % 'full' is not supported yet (have added, but haven't test)
vbhemopt.M = 1; % number of mixture components in each GMM for an emission (should always be 1)

% check validity
if vbhemopt.v0 <= D-1
  error('v0 not large enough...should be > D-1');
end

% convert list of HMMs into an H3M
h3m_b = hmms_to_h3m_hem(hmms, vbhemopt.emit.covar_type, vbhemopt.use_post);


% check seed is set
if isempty(vbhemopt.seed)
  error('vbhemopt.seed has not been set. This is required to obtain reproducible results.');
end


% check for estimating hyps automatically
if iscell(vbhemopt.learn_hyps) || (vbhemopt.learn_hyps == 1)
  do_learn_hyps = 1;
else
  do_learn_hyps = 0;
end

VERBOSE_MODE = vbhemopt.verbose;

%% run for multiple K, multiple S %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the base condition is single K & single S, others just call this base
% condition multiple times

if (length(K)>1) || (length(S) > 1)
  % storage
  out_all = cell(length(K),length(S));
  LLk_all = zeros(length(K),length(S));
  KS_all  = cell(length(K),length(S));
  
  % call learning for each value of K and S
  for ki = 1:length(K)
    for si = 1:length(S)
      myK = K(ki);
      myS = S(si);
      
      vbhemopt1 = vbhemopt;
      vbhemopt1.showplot = 0;
      if (VERBOSE_MODE)
        fprintf('-- VBHEM K=%d, S=%d --\n', myK, myS);
      end
          
      % call learning with a single K and single S
      out_all{ki,si} = vbhem_cluster(hmms, myK, myS, vbhemopt1);
      
      % correct for multiple parameterizations
      LLk_all(ki,si) = out_all{ki,si}.LL + gammaln(myK+1) + gammaln(myS+1);
      % store K,S
      KS_all{ki,si} = [myK, myS];
    end
  end
  
  % get K with max data likelihood
  [maxLLk,indk] = max(LLk_all(:));
  
  % return the best model
  h3m             = out_all{indk};
  h3m.model_LL    = LLk_all;
  h3m.model_KS    = KS_all;
  h3m.model_bestK = KS_all{indk}(1);
  h3m.model_bestS = KS_all{indk}(2);
  h3m.model_all   = out_all;
  
  if VERBOSE_MODE >= 1
    if (vbhemopt.prune_hmms)
      myS = cellfun(@(x) length(x.pdf), h3m.hmms);
      tmpstr = sprintf('; after pruning K=%d, S=%s ', h3m.K, mat2str(myS));
    else
      tmpstr = ' ';
    end
    fprintf('** overall best model: K=%d; S=%d; L=%g%s**\n', ...
      h3m.model_bestK, h3m.model_bestS, maxLLk, tmpstr);
    if strcmp(vbhemopt.initmode, 'auto')
      fprintf('best init was %s',  h3m.Initmodes);
    end
    if (do_learn_hyps)
      fprintf('  estimated hyps: ');
      vbhmm_print_hyps({h3m.learn_hyps.hypinfo.optname}, h3m.learn_hyps.vbopt);
      fprintf('\n');
    end
  end
  
  h3m_out = h3m;

else
  %% run for a single S and single K %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if ~strcmp(vbhemopt.initmode, 'auto')
    % run VBHEM clustering
    
    [h3m_out] = vbhem_h3m_c(h3m_b, vbhemopt);
    
  else
    %% auto initialization - try these methods
    
    if isfield(vbhemopt,'initmodes')
      initmodes = vbhemopt.initmodes;
      initopts = vbhemopt.initopts;
    else
      
      initmodes = { 'baseem', 'gmmNew','wtkmeans'};
      initopts = { 'u', 'r0', 'r0'};
    end
    
    
    for ii=1:length(initmodes)
      tic;
      vbhemopt3 = vbhemopt;
      vbhemopt3.initmode = initmodes{ii};
      vbhemopt3.initopt.mode = initopts{ii};
      fprintf('auto initialization: trying %s: ', vbhemopt3.initmode);
      h3m_out_trials{ii} = vbhem_h3m_c(h3m_b, vbhemopt3);
      h3m_out_trials{ii}.Initmodes = initmodes{ii};
      trials_LL(ii) = h3m_out_trials{ii}.LL;
      h3m_out_trials{ii}.time = toc;
    end
    
    % get best method
    [maxLL, ind] = max(trials_LL);
    h3m_out = h3m_out_trials{ind};
    if vbhemopt.keep_best_random_trial
      h3m_out.h3m_out_trials = h3m_out_trials;
    end
    
    fprintf('best init was %s; LL=%g\n', initmodes{ind}, trials_LL(ind));
    
    % correct for multiple parameterization
    h3m_out.LL_orignal = h3m_out.LL;
  end
  
  % post processing
  % change "hmm" to "hmms"
  h3m_out.hmms = h3m_out.hmm;
  h3m_out = rmfield(h3m_out, 'hmm');
  
  % perform pruning and sorting
  if vbhemopt.prune_hmms
    [h3m_out, emptyinds] = vbh3m_remove_empty(h3m_out, 1.0, VERBOSE_MODE);
  end
  
  % sort ROIs
  if ~isempty(vbhemopt.sortclusters)
    for j = 1:length(h3m_out.hmms)
      tmphmm = h3m_out.hmms{j};
      tmphmm = vbhmm_standardize(tmphmm, vbhemopt.sortclusters);
      h3m_out.hmms{j} = tmphmm;
    end
  end
  
end  %END IF K & S

h3m_r = h3m_out;

% append the options
h3m_r.vbhemopt = vbhemopt;

% append version numbers
h3m_r.emhmm_version = emhmm_version();
h3m_r.matlab_version = version();
  

function hemopt = setdefault(hemopt, field, value)
if ~isfield(hemopt, field)
  hemopt.(field) = value;
end

