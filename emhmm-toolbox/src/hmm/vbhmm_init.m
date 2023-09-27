function [mix_t] = vbhmm_init(datai,K,ini)
% vbhmm_init - initialize vbhmm (internal function)
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-01-13
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong

% 2017-06-14: epsilon changed to row format
% 2017-08-21: changed hyp names to alpha0, etc.
% 2017-11-26: bug fix - handle special cases when N<=K and K=1
% 2018-11-26: bug fix - specify 'Start' as 'randSample' (which was the default previously) to be consistent with previous Matlab versions.
% 2020-03-06: v0.77 - handle more cases when gmdistribution.fit fails
%                    (zero variance in covariance matrix, no such function)
% 2021-02-24: v0.78 - fix_clusters during initialization.

VERBOSE_MODE = ini.verbose;

%data = convert_to_matrix(datai);
data = cat(1,datai{:});
[N,dim] = size(data);

switch(ini.initmode)
  
  %% initialize GMM with random initialization
  case 'random'
    
    if (K==1) 
      % check if just one component, which is a Gaussian
      mix.PComponents = [1];
      mix.mu    = mean(data);
      mix.Sigma = cov(data);
      
    elseif (N<=K)
      % check if enough data to learn a GMM
      % if not, then
      %   select point as the means and iid variance
      Nextra = K-N;
      tmp = [ones(1,N), 0.000001*ones(1,Nextra)];      
      mix.PComponents = tmp / sum(tmp);      
      mix.mu    = [data; zeros(Nextra, dim)];
      tmp = var(data);
      mix.Sigma = mean(tmp)*repmat(eye(dim), [1 1 K]);
      
    else
      % enough data, so run GMM
      try
        warning('off', 'stats:gmdistribution:FailedToConverge');
        
        if ~isempty(ini.random_gmm_opt)
          if (VERBOSE_MODE >= 3)
            fprintf('random gmm init: passing random_gmm_opt = ');
            ini.random_gmm_opt
          end
        end
        
        % use the Matlab fit function to find the GMM components
        % 2017-01-21 - ABC - added ability to pass other options (for Antoine)
        % 2018-11-26 - v0.74 - to be consistent with previous versions, specify 'Start' as 'randSample' (which was the default previously)
        mix = gmdistribution.fit(data,K, ini.random_gmm_opt{:}, ...
          'Options', struct('TolFun', 1e-5), 'Start', 'randSample');
        
      catch ME
        % v0.77 - handle different exceptions better
        switch(ME.identifier)
          case 'stats:gmdistribution:IllCondCovIter'
            if (VERBOSE_MODE >= 2)
              fprintf('using shared covariance');
            end
            % 2018-02-01: added regularize option for stability.
            mix = gmdistribution.fit(data,K, 'SharedCov', true, 'Regularize', 1e-10, ...
              'Options', struct('TolFun', 1e-5), 'Start', 'randSample');
            %mix = gmdistribution.fit(data,K, 'SharedCov', true, 'Options', struct('TolFun', 1e-5));
            
            % SharedCov -> SharedCovariance
          
          case 'stats:gmdistribution:ZeroVariance'
            ME
            fprintf('data has zero variance in one dimension.\ntry using "vbopt.random_gmm_opt = {''Regularize'', 0.0001}".\n')
            error('data has zero variance');
            
          case 'MATLAB:UndefinedFunction'          
            % otherwise use our built-in function (if Stats toolbox is not available)
            if (VERBOSE_MODE >= 2)
              fprintf('using built-in GMM');
              warning('vbhmm_init:nostatstoolbox', 'Stats toolbox is not available -- using our own GMM code. Results may be different than when using the Stats toolbox.');
              warning('off', 'vbhmm_init:nostatstoolbox');
            end
            gmmopt.cvmode = 'full';
            gmmopt.initmode = 'random';
            gmmopt.verbose = 0;
            gmmmix = gmm_learn(data', K, gmmopt);
            mix.PComponents = gmmmix.pi(:)';
            mix.mu = cat(2, gmmmix.mu{:})';
            mix.Sigma = cat(3, gmmmix.cv{:});
            
          otherwise
            ME
            error('unhandled error with gmdistribution.fit')            
        end
      end
    end
    
    %% initialize GMM with given GMM
  case 'initgmm'
    mix.Sigma = cat(3, ini.initgmm.cov{:});    
    mix.mu = cat(1, ini.initgmm.mean{:});    
    if size(mix.mu, 1) ~= K
      error('bad initgmm dimensions -- possibly mean is not a row vector');
    end
    
    mix.PComponents = ini.initgmm.prior(:)';
  
    %% initialize GMM with component splitting
  case 'split'
    gmmopt.cvmode = 'full';
    gmmopt.initmode = 'split';
    gmmopt.verbose = 0;
    gmmmix = gmm_learn(data', K, gmmopt);
        
    mix.PComponents = gmmmix.pi(:)';
    mix.mu = cat(2, gmmmix.mu{:})';
    mix.Sigma = cat(3, gmmmix.cv{:});
    
    %% initialize with HMM
  case 'inithmm'
    % handled later.
    
  otherwise
    error('bad initmode');
end

mix_t = {};
mix_t.dim = dim;
mix_t.K = K;
mix_t.N = N;

% setup hyperparameters
mix_t.alpha0   = ini.alpha0;
mix_t.epsilon0 = ini.epsilon0;
if length(ini.mu0) ~= dim
  error(sprintf('vbopt.mu0 should have dimension D=%d', dim));
end
mix_t.m0       = ini.mu0;
mix_t.beta0    = ini.beta0;
if numel(ini.W0) == 1
  % isotropic W
  mix_t.W0       = ini.W0*eye(dim);  % 2016/04/26: ABC BUG Fix: was inv(ini.W0*eye(dim))
  mix_t.W0mode   = 'iid';
else
  % diagonal W
  if numel(ini.W0) ~= dim
    error(sprintf('vbopt.W0 should have dimension D=%d for diagonal matrix', dim));
  end
  mix_t.W0       = diag(ini.W0);
  mix_t.W0mode   = 'diag';
end
if ini.v0<=dim-1
  error('v0 not large enough');
end
mix_t.v0       = ini.v0; % should be larger than p-1 degrees of freedom (or dimensions)
mix_t.W0inv    = inv(mix_t.W0);

% setup model (M-step)
if strcmp(ini.initmode, 'inithmm')
  % use inithmm  
  mix_t.alpha   = ini.inithmm.varpar.alpha;
  mix_t.epsilon = ini.inithmm.varpar.epsilon;
  mix_t.beta    = ini.inithmm.varpar.beta;
  mix_t.v       = ini.inithmm.varpar.v;  
  mix_t.m       = ini.inithmm.varpar.m;
  mix_t.W       = ini.inithmm.varpar.W;
  
else
  % use GMMs
  mix_t.Nk = N*mix.PComponents';
  mix_t.Nk2 = N*repmat(1/K,K,1);
  
  mix_t.xbar = mix.mu';
  
  if size(mix.Sigma,3) == K
    mix_t.S = mix.Sigma;
  elseif (size(mix.Sigma,3) == 1)
    % handle shared covariance
    mix_t.S = repmat(mix.Sigma, [1 1 K]);
  end
  
  % 2017-01-21: handle diagonal Sigma (for Antoine)
  % for diagonal covarainces, one of the dimensions will have size 1
  if ((size(mix_t.S,1) == 1) || (size(mix_t.S,2)==1)) && (dim > 1)
    oldS = mix_t.S;
    mix_t.S = zeros(dim, dim, K);
    for j=1:K
      mix_t.S(:,:,j) = diag(oldS(:,:,j)); % make the full covariance
    end
  end
  
  mix_t.alpha = mix_t.alpha0 + mix_t.Nk2;
  for k = 1:K
    mix_t.epsilon(k,:) = mix_t.epsilon0 + mix_t.Nk2;
  end
  
  % 2021-02-24: v0.78 ABC: if fixed clusters, then don't update them for initialization
  if ini.fix_clusters
    mix_t.m = mix_t.xbar;   % specified means
    % the covariance will be W / (v-d-1) = S
    mix_t.v = (dim+2)*ones(K,1);
    for j=1:K
      mix_t.W(:,:,j) = inv(mix_t.S(:,:,j));
    end
    mix_t.beta = ones(K,1);
  else
    
    mix_t.beta = mix_t.beta0 + mix_t.Nk;
    mix_t.v = mix_t.v0 + mix_t.Nk + 1;   % 2016/04/26: BUG FIX (add 1)
    mix_t.m = ((mix_t.beta0*mix_t.m0)*ones(1,K) + (ones(dim,1)*mix_t.Nk').*mix_t.xbar)./(ones(dim,1)*mix_t.beta');
    mix_t.W = zeros(dim,dim,K);
    for k = 1:K
      mult1 = mix_t.beta0.*mix_t.Nk(k)/(mix_t.beta0 + mix_t.Nk(k));
      diff3 = mix_t.xbar(:,k) - mix_t.m0;
      mix_t.W(:,:,k) = inv(mix_t.W0inv + mix_t.Nk(k)*mix_t.S(:,:,k) + mult1*diff3*diff3');
    end
  end
end

mix_t.C = zeros(dim,dim,K);
mix_t.const = dim*log(2);
mix_t.const_denominator = (dim*log(2*pi))/2;

%mix_t.maxIter = ini.maxIter;
%mix_t.minDiff = ini.minDiff;

%convert struct to matrix
%function out = convert_to_matrix(data)
%out = []; j = 1;
%for i = 1:size(data,1)
%    tout = data{i,1};
%    k = size(tout,1);
%    out(j:j+k-1,:) = tout;
%    j = j + k;
%end