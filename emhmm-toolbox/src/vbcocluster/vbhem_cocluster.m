function [vbco]= vbhem_cocluster(hmms, K, S, vbhemopt)
% vbhem_cocluster - co-clustering into groups using variational Bayesian HEM
%
% perform co-clustering and select the number of clusters automatically
%
% [vbco]= vbhem_cocluster(hmms, K, S, vbhemopt)
%
% INPUTS
%      hmms = [Stimuli(Trial) x Subject] cell array of HMMs, learned with vbhmm_learn.
%
%         K = number or vecter of cluster (required)  
%         S = number number or vecter of states (ROIs) in one HMM     
%          vector: automatically selects the number of hidden states from the given values. 
%                  The model K,S with highest log-likelihood is selected.
%             
%  vbhembopt = structure containing other options as below:
%
%   VB hyper-parameters: each hmm in h3m_b has the same initial hyper-parameters
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
%     vbhemopt.eta0  = Dirichlet distribution concentration parameter on initial state probabilities (default=1)
%     vbhemopt.epsilon0 = Dirichlet distribution for rows of the transition matrix (default=1). 
%                     The meaning is similar to alpha, but for the transition probabilities.
%     vbhemopt.m0   = prior mean - it should be dimension D;
%                   for D=2 (fixation location), default=[256;192] -- typically it should be at the center of the image.
%                   for D=3 (location & duration), default=[256,192,250]
%     vbhemopt.W0   = size of the inverse Wishart distribution (default=0.005).
%                   If a scalar, W is assumed to be isotropic: W = vbhemopt.W*I.
%                   If a vector, then W is a diagonal matrix:  W = diag(vbhemopt.W);
%                   This determines the covariance of the ROI prior.
%                   For example W=0.005 --> covariance of ROI is 200 --> ROI has standard devation of 14.
%     vbhemopt.v0   = dof of inverse Wishart, v > D-1. (default=5) -- larger values give preference to diagonal covariance matrices.
%     vbhemopt.lambda0  = Wishart concentration parameter -- large value encourages use of
%                   prior mean & W, while small value ignores them. (default=1)
%                   (same as beta0 in vbhmm_learn)
%
%   VBHEM Algorithm parameters
%     vbhemopt.initmode     = initialization method (default='baseem')
%                             follow the initialization ways in VHEM, using
%                             the parameters in the h3m to update the
%                             hyperparameters.
%                            'random' - randomly set assignment varibles z and phi 
%                            'baseem' - randomly select a set of base
%                            emissions
%                            'inith3m' - specify a H3M as the initialization
%                            'gmmNew', 'gmmNew2' - initialize emissions
%                            using HEM-G3M estimated
%                             'gmmNew_r' - combine gmmNew and random,
%                             using the information of mapping in the HEM-G3M
%     vbhemopt.trials    = number of trials for 'random' initialization (default=50)
%                          (more trials will increase the memeory usage. if you have limited memory, 
%                           you can run fewer trials at a time, and change the seed each time).
%     vbhemopt.maxIter      = max number of iterations (default=100)
%     vbhemopt.minDiff      = tolerence for convergence (default=1e-5)
%     vbhemopt.showplot     = show plots (default=0)
%     vbhemopt.sortclusters = '' - no sorting [default]
%                            'e' - sort ROIs by emission probabilites
%                            'p' - sort ROIs by prior probability
%                            'f' - sort ROIs by most-likely fixation path [default]
%                              (see vbhmm_standardize for more options)
%     vbhemopt.verbose     = 0 - no messages
%                          = 1 - a few messages showing progress [default]
%                          = 2 - more messages
%                          = 3 - debugging
%     vbhemopt.keep_K_models = 0 - don't keep models for each K [default]
%                              1 - keep models for each K
%
% OUTPUT
%   vbco                 = structure containing the hmms in the reduced model and related information
%   vbco.K               = number of clusters
%   vbco.label           = labels for each subject
%   vbco.groups          = indices for each group
%   vbco.group_size      = size of each group
%   vbco.L_elbo          = log-likelihood (lower bound) for each stimuli, subject, and cluster
%   vbco.cogroup_hmms{j} = group HMMs structure for each stimuli {ST x Kr}
%                          (this structure compatible with output of VHEM
%                          .omega      = omega [Kr x 1]
%                          .Z          = assignment variables [N x Kr]
%                          .label      = labels for each subject
%                          .groups     = indices for each group
%                          .group_size = size of each group
%                          .hmms       = HMM cell array {1 x Kr}
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2021-07-29
% Antoni B. Chan, Janet H. Hsiao, Lan Hui
% City University of Hong Kong, University of Hong Kong

% 2021-07-29: v0.80 - Hui - initial version
%                   - ABC - optimize storage
%                         - added "keep_K_models"
%                         - add version numbers

if nargin<3
  S = [];
end

if nargin<4  
  vbhemopt = struct;
end

% number of states in a group HMM
sizeall = size(hmms);
stimulisize = sizeall(1);
subjectsize = sizeall(2);

%% set default values
vbhemopt = setdefault(vbhemopt, 'verbose', 1);
VERBOSE_MODE = vbhemopt.verbose;

% number of clusters
vbhemopt.K = K;
% number of states 
% for each stimulus, check for empty HMMs
hmms_validmask = zeros(stimulisize, subjectsize);
hmms_validinds = {};
for i=1:stimulisize
  tmp = ~cellfun(@isempty, hmms(i,:));
  hmms_validmask(i,:) = tmp;
  hmms_validinds{i} = find(tmp);
end

vbhemopt = setdefault(vbhemopt, 'remove_empty', 1);  % remove empty ROIs
% v0.73 - remove empty states
% do this before setting number of states
if (vbhemopt.remove_empty)
    
  fprintf('Checking input HMMS: ');
  for i=1:size(hmms,1)
    for j=1:size(hmms,2)
      if ~isempty(hmms{i,j})
        [hmms{i,j}, zi] = vbhmm_remove_empty(hmms{i,j}, 0, 1e-3);
        if ~isempty(zi)
          fprintf('(%d,%d) removed states', i,j);
          fprintf(' %d', zi);
          fprintf('; ');
        end
      end
    end
  end
  fprintf(' done\n');
end

vbhemopt.S = [];
if isempty(S)
    % iterate through all stimuli (trials)
    for i=1:stimulisize 
      v = hmms_validinds{i};
      % get median number
      hmmsS = cellfun(@(x) length(x.prior), hmms(i,v));
      tmpS = ceil(median(hmmsS)); % if in the middle, then take the larger S

      if (VERBOSE_MODE > 1)
        fprintf('using median number of states: %d\n', tmpS);
      end

      vbhemopt.S  = [vbhemopt.S; tmpS];
    end
else
    vbhemopt.S = S;
end


% assume the data has same dimensions
indvaild = find(hmms_validmask==1,1);
D = size(hmms{indvaild}.pdf{1}.mean,2);
switch(D)    
  case 2
    defmu = [256;192];
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
vbhemopt = setdefault(vbhemopt, 'trials', 50);     % number of trials to run  in vhem is 'numtrials'
vbhemopt = setdefault(vbhemopt, 'termmode', 'L');   % rule to terminate EM iterations
vbhemopt = setdefault(vbhemopt, 'termvalue', 1e-5);
vbhemopt = setdefault(vbhemopt, 'max_iter', 500);    % max number of iterations
vbhemopt = setdefault(vbhemopt, 'min_iter', 1);     % min number of iterations
vbhemopt = setdefault(vbhemopt, 'sortclusters', 'f');
vbhemopt = setdefault(vbhemopt, 'initmode',  'baseem');  % initialization mode
vbhemopt = setdefault(vbhemopt, 'minDiff',   1e-5);
vbhemopt = setdefault(vbhemopt, 'showplot',  0);
vbhemopt = setdefault(vbhemopt, 'random_gmm_opt', {});
vbhemopt = setdefault(vbhemopt, 'keep_ini_hmms', 0);  % remove empty ROIs


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

vbhemopt = setdefault(vbhemopt, 'keep_K_models', 0);

vbhemopt = setdefault(vbhemopt, 'calc_LLderiv', 0);
vbhemopt = setdefault(vbhemopt, 'keep_suboptimal_hmms', 0);
vbhemopt = setdefault(vbhemopt, 'minimizer', 'minimize-cg');
vbhemopt = setdefault(vbhemopt, 'canuseMEX', 1);
vbhemopt = setdefault(vbhemopt, 'verbose_prefix', '');
vbhemopt = setdefault(vbhemopt, 'keep_best_random_trial', 0);
vbhemopt = setdefault(vbhemopt, 'seed', []);
vbhemopt = setdefault(vbhemopt, 'use_post', 1);


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
vbhemopt.M = 1; % number of mixture components in each GMM for an emission (should always be 1)

% check validity
if vbhemopt.v0 <= D-1
  error('v0 not large enough...should be > D-1');
end

% convert list of HMMs into an H3M
%h3m_b = hmms_to_h3m_hem(hmms, 'diag', vbhemopt.use_post);
h3m_bs = [];
for i=1:stimulisize
    H3M = hmms_to_h3m(hmms(i,:), vbhemopt.emit.covar_type);
    h3m_bs = [h3m_bs; H3M];
end

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

vbhemopt.stimulisize = stimulisize;
vbhemopt.subjectsize = subjectsize;



% if strcmp(vbhemopt.initmode,'vbh3m')
%     vbhemopt.hmms = hmms;
% end

%% run for multiple K, single S %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%- consider different conditions for K & S

if (length(K)==1)
  % initialize with vhem coclustering
  if strcmp(vbhemopt.initmode,'vbh3m')
    [Agroup_hmms] = inivb_vhemco(hmms,vbhemopt,K,S);
  end
  
  vbhemopt2 = vbhemopt;
  vbhemopt2.VHintial = Agroup_hmms;
  clear Agroup_hmms
  
  [h3m_out] = vbhem_coh3m_c1(h3m_bs, vbhemopt2, K);
  
else
  
  len_K = length(K);
  Ah3m_out = cell(len_K,1);
  LLk_all = zeros(1,len_K);
  
  % v0.80 - ABC - make into normal for loop
  %               (parfor happens in vhem_cocluster, and vbhem_coh3m_c1)
  for k = 1:len_K
    %  [Ah3m_out{k}] = vbhem_coh3m_c1(h3m_bs, vbhemopt,K(k),S);
    
    if strcmp(vbhemopt.initmode,'vbh3m')
      [Agroup_hmms] = inivb_vhemco(hmms,vbhemopt,K(k),S);
    end
    
    vbhemopt2 = vbhemopt;
    vbhemopt2.VHintial = Agroup_hmms;
    
    %clear Agroup_hmms
    
    [Ah3m_out{k}] = vbhem_coh3m_c1(h3m_bs, vbhemopt2, K(k));
    
    LLk_all(k) = Ah3m_out{k}.LL;
  end
  
  LLk_all = LLk_all  + gammaln(K+1);
  
  
  % get K with max data likelihood
  [maxLLk,indk] = max(LLk_all);
  
  % return the best model
  h3m_out             = Ah3m_out{indk};
  h3m_out.model_LL    = LLk_all;
  h3m_out.model_k     = K;
  h3m_out.model_bestK = K(indk);
  if (vbhemopt.keep_K_models)
    h3m_out.model_all   = Ah3m_out;
  end
  
  
  if VERBOSE_MODE >= 1
    
    fprintf('best model: K=%d; L=%g\n', K(indk), maxLLk);
    if strcmp(vbhemopt.initmode, 'auto')
      fprintf('best init was %s',  h3m_out.Initmodes);
    end
    if (do_learn_hyps)
      fprintf('Hyperparameters:\n');
      vbhmm_print_hyps({h3m_out.learn_hyps.hypinfo.optname}, h3m_out.learn_hyps.vbopt);
      fprintf('\n');
    end
  end
  
end

%% prunning out hmm with small prior
if VERBOSE_MODE >= 1
  fprintf('Pruning model:\n');
end

for t= 1: length(h3m_out.h3m_rs)
  % since group is same, so emptyinds is the same
  [h3m_out.h3m_rs{t}] = vbh3m_remove_empty_new(h3m_out.h3m_rs{t},vbhemopt.sortclusters);
end

h3m_out.groups = h3m_out.h3m_rs{1}.groups;
h3m_out.group_size = h3m_out.h3m_rs{1}.group_size;
h3m_out.label = h3m_out.h3m_rs{1}.label;
h3m_out.K = h3m_out.h3m_rs{1}.K;

for j=1:length(h3m_out.h3m_rs)
  for k=1:h3m_out.K   
    allS(k,j) = length(h3m_out.h3m_rs{j}.hmm{k}.prior);
  end
end

if VERBOSE_MODE >= 1
  fprintf('after pruning: K=%d; S=%s\n',h3m_out.K, mat2str(allS));
end

%h3m_out.L_elbo = h3m_out.h3m_rs{1}.L_elbo;
%h3m_r = h3m_out;
% append the options
%h3m_r.vbhemopt = vbhemopt;
 
%% build the output structure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vbco = h3m_out;


% rename field h3m_rs to cogroup_hmms
vbco.cogroup_hmms = vbco.h3m_rs;
vbco = rmfield(vbco, 'h3m_rs'); 

% change hmm -> hmms
for i=1:length(vbco.cogroup_hmms)
  tmp = vbco.cogroup_hmms{i};
  tmp.hmms = tmp.hmm;
  tmp = rmfield(tmp, 'hmm'); 
  %tmp = rmfield(tmp, 'L_elbo');
  vbco.cogroup_hmms{i} = tmp;
end

% append options / version
vbco.vbhemopt = vbhemopt;
vbco.emhmm_version = emhmm_version();
vbco.matlab_version = version();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hemopt = setdefault(hemopt, field, value)
if ~isfield(hemopt, field)
  hemopt.(field) = value;
end

