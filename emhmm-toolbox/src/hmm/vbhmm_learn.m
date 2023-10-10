% vbhmm_learn - learn HMM with variational Bayesian EM algorithm
%
% [hmm,L] = vbhmm_learn(data,K,vbopt)
%
% INPUTS
%   data = Nx1 cell array, each element is a fixation sequence: 
%             data{i} = [TxD] matrix, where each row is a fixation (x,y) or (x,y,t)
%          N = number of sequences, T = sequence length (can be different for each i)
%          D = dimension: 2 for fixation location (x,y); 3 for location & duration (x,y,d). Duration should be in milliseconds.
%          
%      K = scalar: number of hidden states (clusters)
%          vector: automatically selects the number of hidden states from the given values. 
%                  The model K with highest log-likelihood is selected.
%
%  vbopt = structure containing other options as below:
%
%   VB hyper-parameters:
%     vbopt.alpha0 = Dirichlet distribution concentration parameter -- large value
%                    encourages uniform prior, small value encourages concentrated prior (default=0.1)
%                    Another way to think of it is in terms of "virtual" samples -- A typical way to 
%                    estimate the probability of something is to count the number of samples that it
%                    occurs and then divide by the total number of samples, 
%                    i.e. # times it occurred / # samples.  
%                    The alpha parameter of the Dirichlet adds a virtual sample to this estimate, 
%                    so that the probability estimate is (# times it occurred + alpha) / # samples.  
%                    So for small alpha, then the model will basically just do what the data says.  
%                    For large alpha, it forces all the probabilities to be very similar (i.e., uniform).
%     vbopt.mu0    = prior mean - it should be dimension D;
%                    for D=2 (fixation location), default=[256;192] -- typically it should be at the center of the image.
%                    for D=3 (location & duration), default=[256;192;150]
%     vbopt.W0     = size of the inverse Wishart distribution (default=0.005).
%                    If a scalar, W is assumed to be isotropic: W = vbopt.W*I.
%                    If a vector, then W is a diagonal matrix:  W = diag(vbopt.W);
%                    This determines the covariance of the ROI prior.
%                    For example W=0.005 --> covariance of ROI is 200 --> ROI has standard devation of 14.
%                   
%     vbopt.v0     = dof of inverse Wishart, v > D-1. (default=5) -- larger values give preference to diagonal covariance matrices.
%     vbopt.beta0  = Wishart concentration parameter -- large value encourages use of
%                    prior mean & W, while small value ignores them. (default=1)
%     vbopt.epsilon0 = Dirichlet distribution for rows of the transition matrix (default=0.1). 
%                      The meaning is similar to alpha, but for the transition probabilities.
%
%     vbopt.seed         = seed number for random number generator (required for reproducible results)
%
%     vbopt.learn_hyps = cell array of hyperparameters to estimate automatically, e.g., {'alpha0', 'epsilon0'}.
%                      = 1 indicates estimate all parameters, i.e., {'alpha0', 'epsilon0', 'v0', 'beta0', 'W0', 'mu0'}.
%                      = 0 indicates do not do hyperparameter estimation [default].
%                        Hyperparameters are estimated for each subject (subjects can have different hyperparameters).
%                        Parameters specified in vbopt will be used as the initial set of hyperparameters.                         
%     vbopt.learn_hyps_batch = Estimate hyperparameters on a batch of subjects (all subjects use the same hyperparameters).
%                            = same format as vbopt.learn_hyps.
%
%     vbopt.hyps_max   = structure containing max values for hyperparameters for numerical stability.
%                        values beyond this will be clipped at the max.
%                   .alpha0   = [default=exp(30)]
%                   .epsilon0 = [default=exp(30)]
%                   .v0       = [default=1e4]
%                   .beta0    = [default=exp(30)]
%                   .W0       = [default=exp(30)]
%     vbopt.hyps_min   = structure containing min values for hyperparameters for numerical stability
%                   .alpha0   = [default=exp(-30)]
%                   .epsilon0 = [default=exp(-30)]
%                   .v0       = [default=exp(-20)+dim-1 ]
%                   .beta0    = [default=exp(-30)]
%                   .W0       = [default=exp(-30)]          
%
%   Visualization & Output options:
%     vbopt.showplot     = show plots (default=0)
%     vbopt.bgimage      = background image to show in the fixation plot [default='']
%     vbopt.verbose      = 0 - no messages
%                        = 1 - a few messages showing progress [default]
%                        = 2 - more messages
%                        = 3 - debugging (also dumps a file if problem is found)
%
%   Fixed ROIs:
%     vbopt.fixed_rois = [] - learn the ROIs from the data (default)
%                      = {[x1, y1, width1, height1], ...} - use a fixed set of specified ROIs
%                          a cell array of ROIs specified by mean (x,y) and size (width,height), 
%                          corresponding to 2 standard deviations from the mean.
%                          Only alpha0 and epsilon0 hyperparameters will be estimated if using learn_hyps.
%
%   Miscellaneous HMM options:
%     vbopt.sortclusters = ''  - no sorting
%                          'e' - sort ROIs by emission probabilites
%                          'p' - sort ROIs by prior probability
%                          'f' - sort ROIs by most-likely fixation path [default]
%                              (see vbhmm_standardize for more options)
%     vbopt.groups       = [N x 1] vector: each element is the group index for a sequence.
%                          each group learns a separate transition/prior, and all group share the same ROIs.
%                          default = [], which means no grouping used
%     vbopt.fix_cov      = fix the covariance matrices of the ROIs to the specified matrix.
%                          if specified, the covariance will not be estimated from the data.
%                          The default is [], which means learn the covariance matrices.
%     vbopt.fix_clusters = 1 - keep Gaussian clusters fixed (don't learn the Gaussians)
%                          0 - learn the Gaussians [default]
%
%   EM Algorithm parameters:
%     vbopt.initmode     = initialization method (default='random')
%                            'random' - initialize emissions using GMM with random initialization (see vbopt.numtrials)
%                            'initgmm' - specify a GMM for the emissions (see vbopt.initgmm)
%                            'split' - initialize emissions using GMM estimated with component-splitting
%                            'inithmm' - specify an HMM as the initialization.
%     vbopt.numtrials    = number of trials for 'random' initialization (default=50)
%     vbopt.random_gmm_opt = for 'random' initmode, cell array of options for running "gmdistribution.fit".
%                            The cell array should contain pairs of the option name and value, which are recognized
%                            by "gmdistribution.fit".
%                            For example, {'CovType','diagonal','SharedCov',true,'Regularize', 0.0001}.
%                            This option is helpful if the data is ill-conditioned for the standard GMM to fit.
%                            The default is {}, which does not pass any options.
%     vbopt.initgmm      = initial GMM for 'initgmm':
%                            initgmm.mean{k} = [1 x D]
%                            initgmm.cov{k}  = [D x D]
%                            initgmm.prior   = [1 x K]                        
%     vbopt.inithmm      = initial HMM, estimated with vbhmm_learn
%     vbopt.maxIter      = max number of iterations (default=100)
%     vbopt.minDiff      = tolerence for convergence (default=1e-5)
%
%
% OUTPUT
%   vb_hmm.prior       = prior probababilies [K x 1]
%   vb_hmm.trans       = transition matrix [K x K]: trans(i,j) is the P(x_t=j | x_{t-1} = i)
%   vb_hmm.pdf{j}.mean = ROI mean [1 x D]
%   vb_hmm.pdf{j}.cov  = covariances [D x D]
%   vb_hmm.LL          = log-likelihood of data
%   vb_hmm.gamma       = responsibilities gamma{i}=[KxT] -- probability of each fixation belonging to an ROI.
%   vb_hmm.M           = transition counts [K x K]
%   vb_hmm.N1          = prior counts [K x 1]
%   vb_hmm.N           = cluster sizes [K x 1] (emission counts)
%
%  for using groups
%   vb_hmm.prior{g}    = prior for group g
%   vb_hmm.trans{g}    = transition matrix for group g
%   vb_hmm.M{g}        = transition counts for group g [K x K]
%   vb_hmm.N1{g}       = prior counts for group g [K x 1]
%   vb_hmm.Ng{g}       = cluster sizes of group g
%   vb_hmm.N           = cluster sizes for all groups [K x 1] (emission counts)
%   vb_hmm.group_ids     = group ids
%   vb_hmm.group_inds{g} = indices of data sequences for group g
%   vb_hmm.group_map     = sanitized group IDs (1 to G)
%  
%  internal variational parameters - (internal variables)
%   vb_hmm.varpar.epsilon = accumulated epsilon
%   vb_hmm.varpar.alpha   = accumulated alpha
%   vb_hmm.varpar.beta    = accumulated beta
%   vb_hmm.varpar.v       = accumulated v
%   vb_hmm.varpar.m       = mean
%   vb_hmm.varpar.W       = Wishart
%
%  for model selection:
%   vb_hmm.model_LL    = log-likelihoods for all models tested
%   vb_hmm.model_k     = the K-values used
%   vb_hmm.model_all   = the best model for each K
%   vb_hmm.model_bestK = the K with the best value.
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-01-13
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong

% These are options only used internally (for testing).
%     vbopt.calc_LLderiv = 1 - calculate LL derivatives wrt hyperparameters
%                        = 0 - do not [default]
%     vbopt.useMEX       = 1 - use MEX implementations [default]
%                        = 0 - do not use MEX
%     vbopt.keep_suboptimal_hmms = 0 - do not keep
%                                = 1 - keep the unique set of HMMs after random trials 
%                                      and before hyp optimziation in vb_hmm.suboptimal_hmms{i}
%                                      this is used by vbhmm_learn_batch.
%     vbopt.keep_best_random_trial = 0 - do not keep
%                                    1 - keep (used for debugging)
%     vbopt.keep_K_models          = 0 - do not keep
%                                    1 - keep best model for each K
%     vbopt.minimizer = 'fminunc'        - requires MATLAB optimization toolbox
%                     = 'minimize-lbfgs' - LBFGS optimizer (by Carl Rasmussen)
%                     = 'minimize-bfgs'  - BFGS optimizer (by Carl Rasmussen) [default]
%                     = 'minimize-cg'    - Conjugate gradient (by Carl Rasmussen)
%     vbopt.verbose_prefix = add string in front of all verbose outputs [default='']

% VERSIONS
% 2016-04-21: - initial version from Tim
% 2016-04-25: - speed optimization, documentation
% 2016-04-26: - fix initialization bugs (vbhmm_init): mix_t.W0 and mix_t.v
% 2016-04-27: - initialize with random trials
%             - expose parameters (minDiff, etc)
% 2016-04-29: - bug fix: calculate correct LL lower-bound (vbhmm_em)
%             - fix precision problem for FB algorithm when T is large (vbhmm_em)
%             - bug fix: xi was normalized by row for each sequence n
%                        it should be normalized by matrix for each time t
% 2016-05-05: - model selection over K
% 2016-05-06: - group data - learn separate transition/prior for each group, but share the same ROIs.
% 2016-08-14: - separate function for running FB algorithm (vbhmm_fb)
%             - save variational parameters for output
% 2016-08-15: - initialize with GMM (cell array)
%             - keep clusters fixed, just train the transition/prior parameters (fix_clusters option)
% 2016-12-09: - initialize using splitting-component GMM.
% 2017-01-11: - added verbose option
% 2017-01-17: - update for fixation duration
%             - bug fix: changed vbopt.mean to vbopt.mu
% 2017-01-21: - added for Antoine: constrain covariance to fixed matrix (fix_cov)
%             - added for Antoine: for initialization with random gmm, 
%                     added passable options for gmdistribution.fit (random_gmm_opt)
% 2017-02-13: - code cleanup for vbhmm_em
%             - move lower-bound calculation to vbhmm_em_lb
% 2017-06-09: - added initialization using existing hmm
% 2017-06-11: - added option to calculate derivatives of LL
% 2017-06-14: - code cleanup (epsilon matrix is in row format)
% 2017-07-31: - added hyperparameter estimation (learn_hyps option).
% 2017-08-02: - speed-up vbhmm_fb using a MEX file
% 2017-08-03: - speed-up vbhmm_em and vbhmm_fb
% 2017-08-04: - added option to return all sub-optimal HMMs
% 2017-08-16: - added "learn_hyps_batch" option for batch learning
% 2017-08-18: - added 'bgimage' option for plotting
% 2017-08-21: - changed hyp naming convention to: mu0, W0, beta0, v0, epsilon0
%             - clip hyps at extreme values to improve stability
% 2017-08-25: - set default minimizer to BFGS
% 2017-11-26: bug fix - vbhmm_init: handle special cases when N<=K and K=1
% 2018-02-01: bug fix - vbhmm_em:   handle L=NaN case
% 2018-02-01: bug fix - vbhmm_init: better handle fallback case when using shared covariance with gmm fit.
% 2018-02-09: v0.72 - check invalid setting of v0<D in vbhmm_learn
%             v0.72 - for hyperparam estimation, add 'W0log' and 'W0isqrt' (same as 'W0')
%             v0.72 - added max/min values for W0
%             v0.72 - xxx removed derivative clipping for extreme hyp values
% 2018-02-13: v0.72 - added new option for min/max values of hyperparameters
% 2018-03-07: v0.72 - when selecting among K, also save best model for each K
% 2018-03-10: v0.72 - added support for groups in vbhmm_fb MEX, and hyp learning
% 2018-04-16: v0.72 - improve stability - symmetrize W and C
% 2018-05-15: v0.72 - warn when stats toolbox is not available
% 2018-06-15: v0.74 - bug fix computing cov when ROI has no fixations
% 2018-11-26: v0.74 - remove temporary HMMs not needed in output of vbhmm_learn:
%                       .learn_hyps.hmm_best_random_trial - use an option to include it 
%                       .learn_hyps.inithmm - this is the same as learn_hyps.vbopt.inithmm
%                       .model_all - use an option to keep the best models for each K 
% 2018-11-26: v0.74 - consistency fix: in vbhmm_init, change initialization method of gmdistribution.fit to match previous MATLAB versions
% 2019-01-12: v0.75 - added "seed" option for initializing rng.
% 2019-02-20: v0.75 - verbose prefix strings, add version number
% 2019-06-14: v0.76 - check validity of vbopt.mu0
% 2020-03-06: v0.77 - better error handling when gmdistribution.fit during vbhmm_init
% 2021-06-03: v0.78 - better handle invalid variational parameters.
% 2021-07-21: v0.78 - added fixed_rois option for specifying fixed ROIs
%                   - fix bug keep_suboptimal_hmms for initgmm mode
% 2021-09-02: v0.80 - put matlab version in output hmm

function [hmm,L] = vbhmm_learn(data,K,vbopt)

if nargin<3
  vbopt = struct;
end

%% check old hyp names, and rename if necessary
OLD_NAMES = {{'mean', 'mu0'}, ...
             {'mu', 'mu0'}, ...
             {'alpha', 'alpha0'}, ...
             {'epsilon', 'epsilon0'}, ...
             {'v', 'v0'}, ...
             {'W', 'W0'}, ...
             {'beta', 'beta0'} };

for i=1:length(OLD_NAMES)
  oldhyp = OLD_NAMES{i}{1};
  newhyp = OLD_NAMES{i}{2};
  if isfield(vbopt, oldhyp)
    warning('DEPRECATED: vbopt.%s has been renamed to vbopt.%s', oldhyp, newhyp);
    vbopt.(newhyp) = vbopt.(oldhyp);
    vbopt = rmfield(vbopt, oldhyp);
  end
end

%% set verbose mode
vbopt = setdefault(vbopt, 'verbose', 1);
VERBOSE_MODE = vbopt.verbose;


%% v0.78: handle fixed ROIs
vbopt = setdefault(vbopt, 'fixed_rois', []);  % set default (learnable ROIs)
if ~isempty(vbopt.fixed_rois)
  % check each setting is not overrided
  if isfield(vbopt, 'fix_clusters') && (vbopt.fix_clusters == 0)
    error('fixed ROIs: vbopt.fix_clusters is set to 0, but it should not be set');
  end
  vbopt.fix_clusters = 1;
  
  % ROIs are specified as: {[x1, y1, width1, height1], ...}
  Kf = length(vbopt.fixed_rois);
  if (VERBOSE_MODE >= 1)
    fprintf('using %d fixed ROIs:\n', Kf);
  end
  
  initgmm.mean = {};
  initgmm.cov = {};
  for i=1:Kf
    myx = vbopt.fixed_rois{i}(1);
    myy = vbopt.fixed_rois{i}(2);
    myw = vbopt.fixed_rois{i}(3);
    myh = vbopt.fixed_rois{i}(4);    
    initgmm.mean{end+1} = [myx, myy];
    initgmm.cov{end+1}  = diag((0.25*[myw, myh]).^2);
    if (VERBOSE_MODE >= 1)
      fprintf('  fixed ROI %d: center=(%g,%g), width=%g, height=%g\n', i, myx, myy, myw, myh);
    end
  end
  initgmm.prior = ones(1,Kf) / Kf;
    
  if isfield(vbopt, 'initgmm')
    error('fixed ROIs: vbopt.initgmm is specified, but it should not be set');
  end    
  vbopt.initgmm = initgmm;
  
  if isfield(vbopt, 'initmode') && ~strcmp(vbopt.initmode, 'initgmm')
    error('fixed ROIs: vbopt.initmode is set, but it should not be set');
  end
  vbopt.initmode = 'initgmm';
  
  % check K mismatch
  if length(K) > 1
    error('more than one K is specified with fixed ROIs. K should be set to the number of fixed ROIs');
  end
  if K ~= Kf
    error('K is not set to the number of fixed ROIs');
  end
  
  % set learnable hyperparaameters
  if isfield(vbopt, 'learn_hyps')
    mylearn_hyps = {'alpha0', 'epsilon0'};
    
    if iscell(vbopt.learn_hyps)
      % check if there are any hyps selected that should not be estimated
      for i=1:length(vbopt.learn_hyps)
        if ~any(strcmp(vbopt.learn_hyps{i}, mylearn_hyps))
          error(sprintf('fixed ROIs: cannot set hyp %s as learnable', vbopt.learn_hyps{i}))
        end
      end      
    
    elseif vbopt.learn_hyps == 1
      % for auto learning, only select alpha0 and epsilon0 for estimation
      vbopt.learn_hyps = mylearn_hyps;
    end
  end
 
  % set sortclusters
  if isfield(vbopt, 'sortclusters') && ~isempty(vbopt.sortclusters)
    error('fixed ROIs: vbopt.sortclusters is set, but it is not suggested to set since it will reorder the ROIs');
  end
  vbopt.sortclusters = '';
  
end




%% set default values
vbopt = setdefault(vbopt, 'alpha0', 0.1);  % the Dirichlet concentration parameter

D = size(data{1}, 2);
switch(D)
  case 2
    defmu = [256;192];
  case 3
    defmu = [256;192;150];
  otherwise
    defmu = zeros(D,1);
    %warning(sprintf('no default mu0 for D=%d', D));
end
vbopt = setdefault(vbopt, 'mu0',   defmu); % hyper-parameter for the mean
vbopt = setdefault(vbopt, 'W0',     .005); % the inverse of the variance of the dimensions
vbopt = setdefault(vbopt, 'beta0',  1);
vbopt = setdefault(vbopt, 'v0',     5);
vbopt = setdefault(vbopt, 'epsilon0', 0.1);
  
vbopt = setdefault(vbopt, 'initmode',  'random');
vbopt = setdefault(vbopt, 'numtrials', 50);
vbopt = setdefault(vbopt, 'maxIter',   100); % maximum allowed iterations
vbopt = setdefault(vbopt, 'minDiff',   1e-5);
vbopt = setdefault(vbopt, 'showplot',  0);
vbopt = setdefault(vbopt, 'bgimage', '');
vbopt = setdefault(vbopt, 'sortclusters', 'f');
vbopt = setdefault(vbopt, 'groups', []);
vbopt = setdefault(vbopt, 'fix_clusters', 0);
vbopt = setdefault(vbopt, 'seed', []);

vbopt = setdefault(vbopt, 'random_gmm_opt', {});
vbopt = setdefault(vbopt, 'fix_cov', []);

% minimum/maximum values of hyperparameters, for numerical stability
vbopt = setdefault(vbopt, 'hyps_max', struct);
hyps_max = vbopt.hyps_max;
hyps_max = setdefault(hyps_max, 'alpha0',   1.0686e+13);  % exp(30)
hyps_max = setdefault(hyps_max, 'epsilon0', 1.0686e+13);  % exp(30)
hyps_max = setdefault(hyps_max, 'v0',       1e4);
hyps_max = setdefault(hyps_max, 'beta0',    1.0686e+13);  % exp(30)
hyps_max = setdefault(hyps_max, 'W0',       1.0686e+13);  % exp(30)
vbopt.hyps_max = hyps_max;

vbopt = setdefault(vbopt, 'hyps_min', struct);
hyps_min = vbopt.hyps_min;
hyps_min = setdefault(hyps_min, 'alpha0',    1.0686e-13);  % exp(-30)
hyps_min = setdefault(hyps_min, 'epsilon0', 1.0686e-13);  % exp(-30)
hyps_min = setdefault(hyps_min, 'v0',       2.0612e-09+D-1);  % exp(-20)+dim-1 - should be okay to ~10000 dimensions
hyps_min = setdefault(hyps_min, 'beta0',    1.0686e-13);  % exp(-30)
hyps_min = setdefault(hyps_min, 'W0',       1.0686e-13);  % exp(-30)
vbopt.hyps_min = hyps_min;

vbopt = setdefault(vbopt, 'learn_hyps', 0);
vbopt = setdefault(vbopt, 'learn_hyps_batch', 0);

vbopt = setdefault(vbopt, 'calc_LLderiv', 0);
vbopt = setdefault(vbopt, 'useMEX', 1);
vbopt = setdefault(vbopt, 'keep_suboptimal_hmms', 0);
vbopt = setdefault(vbopt, 'keep_best_random_trial', 0);
vbopt = setdefault(vbopt, 'keep_K_models', 0);
vbopt = setdefault(vbopt, 'minimizer', 'minimize-bfgs');
vbopt = setdefault(vbopt, 'verbose_prefix', '');


% check validity
if vbopt.v0 <= D-1
  error('v0 not large enough...should be > D-1');
end
if ~isempty(vbopt.fix_cov)
  if any(size(vbopt.fix_cov) ~= [D D])
    error('fix_cov covariance is not correct size');
  end
end

% v0.76 - check validity of mu0
tmp0 = size(vbopt.mu0);
if (tmp0(2) ~= 1) && (tmp0(1) ~= D)
  error('vbopt.mu0 should be a D-dim column vector');
end



% v0.75 - check seed is set
if isempty(vbopt.seed)
  error('vbopt.seed has not been set. This is required to obtain reproducible results.');
end

% check K for inithmm
if strcmp(vbopt.initmode, 'inithmm')
  newK = length(vbopt.inithmm.pdf);
  if all(newK ~= K)
    warning('the K passed to vbhmm_learn is different than that in the initial HMM. Using the K from the initial HMM.')
    K = newK;
  end
end

% check for estimating hyps automatically
if iscell(vbopt.learn_hyps) || (vbopt.learn_hyps == 1)
  do_learn_hyps = 1;
else
  do_learn_hyps = 0;
end


%% run for multiple K %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(K)>1
  % turn off plotting
  vbopt2 = vbopt;
  vbopt2.showplot = 0; 
  
  % call learning for each value of K
  out_all = cell(1,length(K));
  LLk_all  = zeros(1,length(K));
  for ki = 1:length(K)
    if (VERBOSE_MODE >= 2)
      fprintf('%s-- K=%d --\n', vbopt.verbose_prefix, K(ki));
    elseif (VERBOSE_MODE == 1)
      fprintf('%s-- vbhmm K=%d: ', vbopt.verbose_prefix, K(ki));
    end
    
    % for initgmm, select one
    if strcmp(vbopt2.initmode, 'initgmm')
      vbopt2.initgmm = vbopt.initgmm{ki};
    end
    
    % call learning with a single k
    out_all{ki} = vbhmm_learn(data, K(ki), vbopt2);
    LLk_all(ki) = out_all{ki}.LL;
  end

  % correct for multiple parameterizations
  % v0.80 - BUG Fix if K is not a row vector
  LLk_all = LLk_all + gammaln(K(:)'+1);
  
  % get K with max data likelihood
  [maxLLk,ind] = max(LLk_all);
  
  % return the best model
  hmm       = out_all{ind};
  hmm.model_LL    = LLk_all;
  hmm.model_k     = K;
  hmm.model_bestK = K(ind);
  % 2018-11-26: v0.74 - option to keep best model for each K
  if (vbopt.keep_K_models)
    hmm.model_all   = out_all;
  end
  L = maxLLk;
  
  if VERBOSE_MODE >= 1
    fprintf('%sbest model overall: K=%d; L=%g\n', vbopt.verbose_prefix, K(ind), maxLLk);
    
    if (do_learn_hyps)
      fprintf('%s  ', vbopt.verbose_prefix);
      vbhmm_print_hyps({hmm.learn_hyps.hypinfo.optname}, hmm.learn_hyps.vbopt);
      fprintf('\n');
    end    
  end
  
  if (vbopt.keep_suboptimal_hmms)
    % append suboptimal HMMs from all Ks
    tmp = {};
    for k=1:length(out_all)
      tmp = [tmp, out_all{k}.suboptimal_hmms];
    end
    hmm.suboptimal_hmms = tmp;
  end
  
  if 0 && vbopt.showplot
    figure
    hold on
    plot(K, LLk_all, 'b.-')
    plot([min(K), max(K)], [maxLLk, maxLLk], 'k--');
    plot(K(ind), maxLLk, 'bo');
    hold off
    grid on
    xlabel('K');
    ylabel('data log likelihood');
  end


else 
  %% run for a single K %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % v0.75 - set the seed here
  if ~isempty(vbopt.seed)
    if (VERBOSE_MODE > 2)
      fprintf('%s+ setting random seed to %d\n', vbopt.verbose_prefix, vbopt.seed);
    elseif (VERBOSE_MODE == 1)
      fprintf('(seed=%d)', vbopt.seed);
    end
    rng(vbopt.seed, 'twister');
  end
  
  switch(vbopt.initmode)
    %%% RANDOM initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'random'
      numits  = vbopt.numtrials;
      
      % if only one ROI, then only one trial is needed since there is only one solution.
      if (K == 1)
        numits = 1;
      end
      
      vb_hmms = cell(1,numits);
      LLall   = zeros(1,numits);
      % run several iterations
      for it = 1:numits                     
        vb_hmms{it} = vbhmm_em(data, K, vbopt);
        LLall(it) = vb_hmms{it}.LL;
        
        if (VERBOSE_MODE == 1)
          fprintf('.');
          if mod(it,10)==0
            fprintf('%d', it);
          end
        elseif (VERBOSE_MODE > 1)
          fprintf('%sTrial %d: LL=%g\n', vbopt.verbose_prefix, it, LLall(it));
        end 
        
        % debugging
        %  vbhmm_plot(vb_hmms{it}, data);             
      end
      
      show_LL_plot = 0;
      % get unique solutions
      if (vbopt.keep_suboptimal_hmms || do_learn_hyps)
        % if percentage change is larger than twice convergence criteria.
        % 2*   because two models could be on either side of the true model.
        % 10*  since convergence might not happen completely
        diffthresh = 2*vbopt.minDiff*10;
      
        [unique_hmm_it, VLL, legs] = uniqueLL(LLall, diffthresh, show_LL_plot);
        
        unique_hmm    = {vb_hmms{unique_hmm_it}};
        unique_hmm_LL = LLall(unique_hmm_it);
      end
      
      
      %% learn hyps after random initializations
      if (do_learn_hyps)
        show_hyp_plot = 0;
        
        % save old values
        LLall_old   = LLall;
        vb_hmms_old = vb_hmms;
        
        % invalidate
        LLall = nan*LLall_old;
        vb_hmms = cell(size(vb_hmms));
        
        for q=1:length(unique_hmm)
          my_it  = unique_hmm_it(q);
          my_hmm = unique_hmm{q};
          my_LL  = unique_hmm_LL(q);
          
          if (vbopt.verbose >= 1)
            fprintf('\n%s(K=%d) optimizing trial %d: %g', vbopt.verbose_prefix, K, my_it, my_LL);
          end
          
          % estimate hyperparameters
          vbopt2 = vbopt;
          vbopt2.initmode = 'inithmm';
          vbopt2.inithmm  = my_hmm;
          vb_hmms{my_it} = vbhmm_em_hyp(data, K, vbopt2);
          LLall(my_it)  = vb_hmms{my_it}.LL;
            
          % visualization
          if (show_LL_plot)
            figure(VLL)
            legs4 = plot(my_it, LLall(my_it), 'rx');
            legend([legs, legs4], {'original LL', 'selected trial', 'optimized LL', 'ignore bounds'});
            drawnow
          end
        end
        
        if (show_hyp_plot)
         % show models
         faceimg = '../face.jpg';
         figure         
         NN = length(unique_hmm);
         for i=1:NN
           myit = unique_hmm_it(i);
           subplot(2,NN,i)
           vbhmm_plot_compact(unique_hmm{i}, faceimg);
           title(sprintf('trial %d\norig LL=%g', myit, unique_hmm{i}.LL));
           subplot(2,NN,i+NN)
           vbhmm_plot_compact(vb_hmms{myit}, faceimg);
           title(sprintf('opt LL=%g', vb_hmms{myit}.LL));
         end
         drawnow
        end
                
         %LLall
      end % end do_learn_hyps
      
      %% choose the best
      [maxLL,maxind] = max(LLall);
      if (VERBOSE_MODE >= 1)        
        fprintf('\n%s(K=%d) best run=%d; LL=%g\n', vbopt.verbose_prefix, K, maxind, maxLL);
      end
      
      
      %LLall
      %maxLL
      %maxind
      hmm = vb_hmms{maxind};
      L = hmm.LL;           
      
      if (do_learn_hyps)
         % check degenerate models and save (DEBUGGING)
         if any(abs(LLall_old./LLall)>10)
           if (VERBOSE_MODE >= 3)
             foo=tempname('.');
             warning('degenerate models detected: saving to %s', foo);
             save([foo '.mat'])
           else
             warning('degenerate models detected');
           end
         end
         
        
        % output final params
        if (VERBOSE_MODE >= 1)
          fprintf('%s(K=%d)  ', vbopt.verbose_prefix, K);
          vbhmm_print_hyps({hmm.learn_hyps.hypinfo.optname}, hmm.learn_hyps.vbopt);
          fprintf('\n');
        end
        
        % 2018-11-26: v0.74 - use option to include best random init.
        % choose the best among the random initialization
        if (vbopt.keep_best_random_trial)
          [maxLL_old,maxind_old] = max(LLall_old);
          tmphmm = vb_hmms_old{maxind_old};
          if ~isempty(vbopt.sortclusters)
            tmphmm = vbhmm_standardize(tmphmm, vbopt.sortclusters);
          end
          hmm.learn_hyps.hmm_best_random_trial = tmphmm;
        end
        
      end
      
      if (vbopt.keep_suboptimal_hmms)
        hmm.suboptimal_hmms = unique_hmm;
      end
      
    %%% Initialize with learned GMM %%%%%%%%%%%%%%%%%%%%%%
    %%% Initialize with component-splitting GMM %%%%%%%%%%
    %%% Initialize with learned HMM %%%%%%%%%%%%%%%%%%%%%
    case {'initgmm', 'split', 'inithmm'}
      
      if ~(do_learn_hyps)
        % do not learn hyps: just run EM
        hmm = vbhmm_em(data, K, vbopt);
      
      else
        % learn hyps...
        % first get the initial hmm, if not specified
        if ~strcmp(vbopt.initmode, 'inithmm')
          hmm_init = vbhmm_em(data, K, vbopt);          
          vbopt2 = vbopt;
          vbopt2.inithmm = hmm_init;
          vbopt2.initmode = 'inithmm';
        else
          % inithmm already specified
          vbopt2 = vbopt;
        end
        % now learn hyps with init
        hmm = vbhmm_em_hyp(data, K, vbopt2);
      end
      
      L = hmm.LL;
      
      % v0.78: keep suboptimal hmms
      if (vbopt.keep_suboptimal_hmms)
        hmm.suboptimal_hmms = {hmm};
      end
      
  end


end

% post processing
% - reorder based on cluster size
if ~isempty(vbopt.sortclusters)
  hmm_old = hmm;
  hmm = vbhmm_standardize(hmm, vbopt.sortclusters);
  %[wsort, wi] = sort(hmm.N, 1, 'descend');
  %hmm = vbhmm_permute(hmm, wi);
end

% plotting
if (vbopt.showplot)
  vbhmm_plot(hmm, data, vbopt.bgimage);
  drawnow
end

% append the options
hmm.vbopt = vbopt;

% append version number
hmm.emhmm_version = emhmm_version();
hmm.matlab_version = version();
  
function vbopt = setdefault(vbopt, field, value)
if ~isfield(vbopt, field)
  vbopt.(field) = value;
end
