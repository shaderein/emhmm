function h3m = initialize_hem_h3m_c(h3m_b,mopt)
% 'r'       random initilization (do not use)
% 'base'    randomly select a subset of the input hmms
% 'highp'   select the input hmms with largest weight
% 'gmmNew'  estimate GMM from pooled emissions - each HMM component is
%           one Gaussian from the GMM.
% 'gmm'     use HEM-GMM on the emission GMMs (pooled together), than
% initialize each HMM cluster center by randomizing the weight of the
% Gaussian, and with random values for the dynamics
% 'gmm_Ad'  similar to 'gmm', but set a strong diagonal component to
% transition matrices
% 'gmm_L2R' similar to 'gmm', but for left to right HMMs
%
% 'baseem' randomly select a set of base emissions
% ---
% H3M Toolbox 
% Copyright 2012, Emanuele Coviello, CALab, UCSD

% ABC 2017-05-18 - fixed 'gmm' mode to better handle H3Ms w/o the same number of states.
%                - some bug fixes
% ABC 2017-05-30 - fixed 'gmmNew' mode to better handle H3Ms w/o the same number of states.
% ABC 2017-08-17 - added full covariance matrices
% 2018-11-23: ABC v0.74 - handle empty HMMs (baseem)
%                 v0.74 - NOTE: still need to do gmmNew and gmmNew2

% mopt.initopt.mode   A-prior initialization
%                     = 'u' - uniform prior and A
%                     = 'r' - random prior and A%
%                     = 'm' - use mapping of ROIs from base to reduced
%                     GMM-initialization
%                     = '0' - use steady-state state probabilities as component weights
%                     = '3' - use average state probability of first 3 time-steps
%                     = 'x' - use uniform weights

h3m_K = h3m_b.K;            

% preprocess data
switch(mopt.initmode)
  % when selecting base HMM components, ensure correct number of states
  case {'base'}
    %%%Remove HMMs w/o correct number of states
    b = 1;
    h3m_b_t = {};
    for m = 1:h3m_K
      aaa = size(h3m_b.hmm{m}.prior,1);
      if aaa == mopt.N
        h3m_b_t.hmm{1,b} = h3m_b.hmm{1,m};
        b = b+1;
      end
    end
    b = b - 1;
    h3m_b_t.K = b;
    omegas = ones(1,h3m_b_t.K);
    h3m_b_t.omega = omegas;
    h3m_b_t.omega = h3m_b_t.omega/sum(h3m_b_t.omega);
    h3m_b = h3m_b_t;
    
  otherwise
    % do nothing
end

mopt.Nv = 1000 * h3m_b.K;

T = mopt.tau;

if isfield(h3m_b.hmm{1}.emit{1},'nin')
    dim = h3m_b.hmm{1}.emit{1}.nin;
else
    dim = size(h3m_b.hmm{1}.emit{1}.centres,2);
end

Kb = h3m_b.K;
Kr = mopt.K;
N = mopt.N;
Nv = mopt.Nv;
M = mopt.M;

switch mopt.initmode
  
    % random (don't use)
  case 'r'
    if strcmp(h3m_b.hmm{1}.emit{1}.covar_type, 'full')
      error('full not supported');
    end
    
    h3m.K = Kr;
    for j = 1 : Kr
      prior = rand(N,1);
      prior = prior/sum(prior);
      A     = rand(N);
      A     = A ./ repmat(sum(A,2),1,N);
      emit = cell(1,N);
      for n = 1 : N
        emit{n} = gmm(dim, M, mopt.emit.covar_type);
        emit{n}.covars = emit{n}.covars .* 100;
        % and modify priors
        emit{n}.priors = rand(size(emit{n}.priors));
        emit{n}.priors = emit{n}.priors / sum(emit{n}.priors);
      end
      
      h3m.hmm{j}.prior = prior;
      h3m.hmm{j}.A     = A;
      h3m.hmm{j}.emit  = emit;
    end
    omega = rand(1,Kr);
    omega = omega/sum(omega);
    
    h3m.omega = omega;
    h3m.LogL = -inf;
    
  case 'baseem'
    h3m.K = Kr;
    
    for j=1:Kr
      % 2018-11-23: v0.74 - skip hmms w/ 0 weight
      %for n=1:N
      n=1;
      while(n<=N)
        randomb = randi(Kb);
        % 2018-11-23: v0.74 - skip hmms w/ 0 weight
        if (h3m_b.omega(randomb) == 0)
          continue
        end
        randomg = randi(length(h3m_b.hmm{randomb}.emit));
        %fprintf('%d,%d; ', randomb, randomg);
        h3m.hmm{j}.emit{n} = h3m_b.hmm{randomb}.emit{randomg};
        n = n+1;
      end
      % old code used uniform
      %h3m.hmm{j}.prior = ones(N,1)/N;
      %h3m.hmm{j}.A     = ones(N,N)/N;
      %[h3m.hmm{j}.prior, h3m.hmm{j}.A] = makeAprior(N, 'u');

      [h3m.hmm{j}.prior, h3m.hmm{j}.A] = makeAprior(N, mopt.initopt.mode);
    end
    
    omega = ones(1,Kr)/Kr;
    h3m.omega = omega;
    h3m.LogL = -inf;
    
    
  case 'base'
    indexes = randperm(Kb);
    h3m.K = Kr;
    for j = 1 : Kr
      
      h3m.hmm{j} = h3m_b.hmm{indexes(j)};
    end
    %omega = rand(1,Kr);
    omega = ones(1,Kr);
    omega = omega/sum(omega);
    h3m.omega = omega;
    h3m.LogL = -inf;
    
    % 2017-08-09 : ABC - k-means initialization [INCOMPLETE]
  case 'kmeans'
    error('kmeans init not fully tested');
    if strcmp(h3m_b.hmm{1}.emit{1}.covar_type, 'full')
      error('full not supported');
    end
    
    % cluster the most-likely center-trajectories
    Tlength = 3;
    CT = zeros(dim, Tlength, Kb);  % center trajectories
    stateT = zeros(Tlength, Kb);     % state trajectory
    emitT  = zeros(Tlength, Kb);    % emission gmm trajectory
    for i=1:Kb
      for t=1:Tlength
        % get state probability 
        if t==1
          px = h3m_b.hmm{i}.prior';
        else
          px = px*h3m_b.hmm{i}.A;
        end
        % select most likely center
        [~,j] = max(px);
        [~,k] = max(h3m_b.hmm{i}.emit{j}.priors);
        CT(:,t,i) = h3m_b.hmm{i}.emit{j}.centres(k,:);
        stateT(t,i) = j;
        emitT(t,i) = k;
      end
    end
    
    % normalize by the first fixation
    CTn = bsxfun(@minus, CT, CT(:,1,:));
    
    CT2 = reshape(CTn,dim*Tlength,Kb)';
    CT2
    [inds,Cs] = kmeans(CT2, Kr); %, 'distance', 'cosine');
    
    % get group members
    for k=1:Kr
      group{k} = find(inds(:)'==k);
    end
    
    showplots = 1;
    
    if (showplots)      
      % plot trajectories
      figure
      colors = 'rgbcymk';
      hold on
      Cstmp = reshape(Cs',[dim, Tlength, Kr]);
      for k=1:Kr
        for iii=group{k}
          plot(CTn(1,:,iii), CTn(2,:,iii), [colors(k) 'x-']);
        end
        plot(Cstmp(1,:,k), Cstmp(2,:,k), [colors(k) 'o-']);
      end
      hold off
      
      % plot clusters
      h3m_b.Z = [];
      h3m_b.LogLs = [];
      h3m_b.LogL = [];
      ghmms = h3m_to_hmms(h3m_b, 'diag');
      for k=1:Kr
        figure
        px = 5;
        py = ceil(length(group{k})/5);
        for i=1:length(group{k})
          subplot(py,px,i)
          vbhmm_plot_compact(ghmms.hmms{group{k}(i)}, '../test_jov/ave_face120.png');
        end
      end
      inds
    end
    
    % [TODO]
    % make h3m    
    h3m.K = Kr;    
    for j=1:Kr      
      h3m.hmm{j}.prior = ones(N,1)/N;
      h3m.hmm{j}.A     = ones(N,N)/N;
      h3m.hmm{j}.emit{i}.type     = 'gmm';  % a GMM with one component
      h3m.hmm{j}.emit{i}.nin      = nin;
      h3m.hmm{j}.emit{i}.ncentres = 1;
      h3m.hmm{j}.emit{i}.priors   = 1;
      h3m.hmm{j}.emit{i}.centres  = hmms{j}.pdf{i}.mean;
      h3m.hmm{j}.emit{i}.covar_type = 'diag';
      h3m.hmm{j}.emit{i}.covars = diag(hmms{j}.pdf{i}.cov)';        
    end
    h3m.omega = ones(1,Kr)/Kr;;
    h3m.LogL = -inf;
    
    
  case 'trick'
    indexes = 1:(Kb/Kr):Kb;
    h3m.K = Kr;
    for j = 1 : Kr
      
      h3m.hmm{j} = h3m_b.hmm{indexes(j)};
    end
    omega = rand(1,Kr);
    omega = omega/sum(omega);
    h3m.omega = omega;
    h3m.LogL = -inf;
    
  case 'highp'
    [foo indexes] = sort(h3m_b.omega,'descend');
    h3m.K = Kr;
    for j = 1 : Kr
      
      h3m.hmm{j} = h3m_b.hmm{indexes(j)};
    end
    omega = ones(1,Kr)/Kr;
    omega = omega/sum(omega);
    h3m.omega = omega;
    h3m.LogL = -inf;

    % use HEM-G3M - for HMMs with single Gaussian components
    case {'gmm2'}
        
      
    % use hem-g3m    
  case {'gmmNew', 'gmmNew2'}
    % gmmNew - estimate N Gaussians for N states - share Gaussians across groups
    % gmmNew2 - estimate N*Kr Gaussians - each group has different Gaussians.
    
    virtualSamples = Nv * Kb;
    iterations = mopt.initopt.iter;
    trials = mopt.initopt.trials;
    % fit a gaussian mixture to the data...
    
    %% OLD CODE - this doesn't handle individual HMMs with different number of states
    if 0
      % if there is too much input, use only some data
      gmms_all = cell(1,Kb*N);
      alpha = zeros(1,Kb*N);
      
      for i = 1 : h3m_b.K
        
        gmms_all( (i-1)*N+1 : (i)*N ) = h3m_b.hmm{i}.emit;
        % alpha( (i-1)*N+1 : (i)*N ) = h3m_b.hmm{i}.prior;
        p = h3m_b.hmm{i}.prior';
        A = h3m_b.hmm{i}.A;
        for t = 1 : 50
          p = p * A;
        end
        alpha( (i-1)*N+1 : (i)*N ) = p;
      end
    end
    
    % 2017-05-30: ABC
    %% NEW CODE: handles different number of ROIs
    getN = @(x) length(x.prior);
    numbc = sum(cellfun(getN, h3m_b.hmm));
    
    gmms_all = cell(1,numbc);
    alpha = zeros(1,numbc);
    gmms_hmm_index = zeros(1,numbc);
    gmms_hmm_emit_index = zeros(1,numbc);
    
    curj = 1;
    for i = 1 : h3m_b.K
      myS = length(h3m_b.hmm{i}.prior);
      newind = curj:(curj+myS-1);
      gmms_all( curj:(curj+myS-1) ) = h3m_b.hmm{i}.emit;
      %p = h3m_b.hmm{i}.prior';
      %A = h3m_b.hmm{i}.A;
      %for t = 1 : 50
      %  p = p * A;
      %end
      % original method
      %p = makeGMMweights(h3m_b.hmm{i}.prior', h3m_b.hmm{i}.A, '0'); % steady-state
      
      p = makeGMMweights(h3m_b.hmm{i}.prior', h3m_b.hmm{i}.A, mopt.initopt.mode);
      
      gmms_hmm_index(newind) = i;
      gmms_hmm_emit_index(newind) = 1:length(newind);
      
      alpha( newind) = p;
      curj = curj+myS;
    end
    
    %%%
    
    alpha = alpha/sum(alpha);
    gmms.alpha = alpha;
    gmms.mix   = gmms_all;
    
    trainMix(length(gmms_all)) = struct(gmms_all{1});
    for pippo = 1 : length(gmms_all)
      trainMix(pippo).mix = gmms_all{pippo};
      
    end
    
    if strcmp(mopt.initmode, 'gmmNew2')
      tmpK = N*Kr; % separate ROIs for each group
      useGauss = reshape(randperm(tmpK),Kr,N);  % [Kr x N] index into ROIs
    else
      tmpK = N; % share ROIs for each group
      useGauss = repmat(1:N,Kr,1);              % [Kr x N] index into ROIs
    end
    
    % ABC 2017-05-18: BUG fix: should be N (usually Kr>N for tagging)
    [reduced_out, post_out, lp_out] = GMM_MixHierEM(trainMix, tmpK, virtualSamples, iterations);
    %[reduced_out] = GMM_MixHierEM(trainMix, Kr, virtualSamples, iterations);
    
    % map from HMM emission to group ROI
    emit_map_ll = cell(1,Kb);
    for i=1:numbc
      emit_map_ll{gmms_hmm_index(i)}(:,gmms_hmm_emit_index(i)) = lp_out(:,i);
    end
        
    % assign base HMMS to reduced HMM using emissions
    hmm_map = zeros(1,Kb);
    hmm_map_emit = cell(1,Kb);
    for i=1:Kb
      % compute sum of LL between group ROIs and HMM emissions
      emit_LL = zeros(1,Kr);
      for j=1:Kr
        emit_LL(j) = sum(sum(emit_map_ll{i}(useGauss(j,:),:)));
      end
      % pick largest sum LL as the group
      [~, myj] = max(emit_LL);
      % map the base ROIs to group ROIs.
      [~, myii] = max(emit_map_ll{i}(useGauss(myj,:),:),[], 1);
      hmm_map(i) = myj;
      hmm_map_emit{i} = myii;                       
    end
    
    % build new transition matrices
    newA = cell(1,Kr);
    newprior = cell(1,Kr);
    newN = zeros(1,Kr);   
    for j=1:Kr
      newA{j} = zeros(N,N);
      newprior{j} = zeros(N,1);
      newN(j) = 0;
    end
    for i=1:Kb      
      mymap = hmm_map_emit{i}; % get the ROI map from base to reduced
      myj = hmm_map(i);      % get the corresponding group
      % for each base hidden state
      for k=1:length(mymap)
        nk = mymap(k);
        % add probability to corresponding reduced hidden state
        newprior{myj}(nk) = newprior{myj}(nk) + h3m_b.hmm{i}.prior(k);
        for l=1:length(mymap)
          nl = mymap(l);
          newA{myj}(nk,nl) = newA{myj}(nk,nl) + h3m_b.hmm{i}.A(k,l);
        end
      end
      newN(myj) = newN(myj)+1;
    end
    % normalize
    for j=1:Kr
      newprior{j} = newprior{j} / sum(newprior{j});
      newA{j} = bsxfun(@rdivide, newA{j}, sum(newA{j},2));
    end
      
    
    
    if 0
      % show a plot
      figure
      subplot(2,1,1)
      axis ij
      hold on
      cols = {'r-', 'b-', 'g-', 'c-', 'y-', 'k-', 'm-', ...
        'r--', 'b--', 'g--', 'c--', 'y--', 'k--', 'm--'};
      for i=1:length(trainMix)
        [~,x] = max(post_out(:,i));
        switch(trainMix(i).mix.covar_type)
          case 'diag'
            tmp = diag(trainMix(i).mix.covars);
          case 'full'
            tmp =trainMix(i).mix.covars;
        end
        plot2D(trainMix(i).mix.centres, tmp, cols{x});
      end
      hold off
      subplot(2,1,2)
      axis ij
      hold on
      for i=1:reduced_out.ncentres
        
        switch(trainMix(i).mix.covar_type)
          case 'diag'
            tmp = diag(reduced_out.covars(i,:));
          case 'full'
            tmp = reduced_out.covars(:,:,i);
        end
        plot2D(reduced_out.centres(i,:), tmp, cols{i});
      end
      hold off
    end
    
    
    h3m.K = Kr;
    for j = 1 : Kr
      
      emit = cell(1,N);
      for n = 1 : N
        myn = useGauss(j, n);        
        reduced_out_use.centres = reduced_out.centres(myn,:);
        switch(h3m_b.hmm{1}.emit{1}.covar_type)
          case 'diag'
            reduced_out_use.covars = reduced_out.covars(myn,:);
          case 'full'
            reduced_out_use.covars = reduced_out.covars(:,:,myn);
        end
        reduced_out_use.ncentres = 1;
        reduced_out_use.priors = 1;
        reduced_out_use.nin = reduced_out.nin;
        reduced_out_use.covar_type = reduced_out.covar_type;
        reduced_out_use.type = reduced_out.type;
        
        emit{n} = reduced_out_use;
        % and randomize priors ...
        emit{n}.priors = rand(size(emit{n}.priors));
        emit{n}.priors = emit{n}.priors / sum(emit{n}.priors);
      end
      
      % ABC - 2017-05-18: Make it random prior/A      
      % original case was random      
      if any(mopt.initopt.mode == 'm')
        h3m.hmm{j}.prior = newprior{j};
        h3m.hmm{j}.A     = newA{j};
      else
        [prior, A] = makeAprior(N, mopt.initopt.mode);
        h3m.hmm{j}.prior = prior;
        h3m.hmm{j}.A     = A;
      end  
      
      
      h3m.hmm{j}.emit  = emit;
    end
    omega = rand(1,Kr);
    omega = omega/sum(omega);
    h3m.omega = omega;
    h3m.LogL = -inf;    
    
  case {'g3k' 'g3n' 'gmm' 'gmm_Ad' 'gmm_L2R' 'g3k_Ad' 'g3n_Ad' 'gmm_A' 'gmm_Au'}
    if strcmp(h3m_b.hmm{1}.emit{1}.covar_type, 'full')
      error('full not supported');
    end
    
    virtualSamples = Nv * Kb;
    iterations = mopt.initopt.iter;
    trials = mopt.initopt.trials;
    % fit a gaussian mixture to the data...
    
    % if there is too much input, use only some data
    
    
    %%% OLD CODE: this doesn't handle individual HMMs with different number of states
    if 0
      gmms_all = cell(1,Kb*N);
      alpha = zeros(1,Kb*N);
      
      for i = 1 : h3m_b.K
        
        gmms_all( (i-1)*N+1 : (i)*N ) = h3m_b.hmm{i}.emit;
        % alpha( (i-1)*N+1 : (i)*N ) = h3m_b.hmm{i}.prior;
        p = h3m_b.hmm{i}.prior';
        A = h3m_b.hmm{i}.A;
        for t = 1 : 50
          p = p * A;
        end
        alpha( (i-1)*N+1 : (i)*N ) = p;
      end
    end
    
    % 2017-05-18: ABC
    %%% NEW CODE: handles different number of ROIs
    getN = @(x) length(x.prior);
    numbc = sum(cellfun(getN, h3m_b.hmm));
    
    gmms_all = cell(1,numbc);
    alpha = zeros(1,numbc);
    curj = 1;
    for i = 1 : h3m_b.K
      myS = length(h3m_b.hmm{i}.prior);
      newind = curj:(curj+myS-1);
      gmms_all( curj:(curj+myS-1) ) = h3m_b.hmm{i}.emit;
      p = h3m_b.hmm{i}.prior';
      A = h3m_b.hmm{i}.A;
      for t = 1 : 50
        p = p * A;
      end
      alpha( newind) = p;
      curj = curj+myS;
    end
    
    %%%
    
    alpha = alpha/sum(alpha);
    gmms.alpha = alpha;
    gmms.mix   = gmms_all;
    
    trainMix(length(gmms_all)) = struct(gmms_all{1});
    for pippo = 1 : length(gmms_all)
      trainMix(pippo).mix = gmms_all{pippo};
      
    end
    
    % only one Gaussian is estimated
    [reduced_out] = GMM_MixHierEM(trainMix, M, virtualSamples, iterations);
    
    switch mopt.initmode
      case {'gmm'}
        
        h3m.K = Kr;
        for j = 1 : Kr
          prior = rand(N,1);
          prior = prior/sum(prior);
          A     = rand(N);
          A     = A ./ repmat(sum(A,2),1,N);
          emit = cell(1,N);
          for n = 1 : N
            emit{n} = reduced_out;
            % and randomize priors ...
            emit{n}.priors = rand(size(emit{n}.priors));
            emit{n}.priors = emit{n}.priors / sum(emit{n}.priors);
          end
          
          h3m.hmm{j}.prior = prior;
          h3m.hmm{j}.A     = A;
          h3m.hmm{j}.emit  = emit;
        end
        omega = rand(1,Kr);
        omega = omega/sum(omega);
        h3m.omega = omega;
        h3m.LogL = -inf;
        
        
      case {'gmm_Au'}
        
        h3m.K = Kr;
        for j = 1 : Kr
          prior = ones(N,1);
          prior = prior/sum(prior);
          
          A = ones(N);
          A     = A ./ repmat(sum(A,2),1,N);
          emit = cell(1,N);
          for n = 1 : N
            emit{n} = reduced_out;
            % and randomize priors ...
            emit{n}.priors = rand(size(emit{n}.priors));
            emit{n}.priors = emit{n}.priors / sum(emit{n}.priors);
          end
          
          h3m.hmm{j}.prior = prior;
          h3m.hmm{j}.A     = A;
          h3m.hmm{j}.emit  = emit;
        end
        omega = rand(1,Kr);
        omega = omega/sum(omega);
        h3m.omega = omega;
        h3m.LogL = -inf;
        
      case {'gmm_Ad'}
        
        h3m.K = Kr;
        for j = 1 : Kr
          prior = rand(N,1);
          prior = prior/sum(prior);
          
          % pick one A form the examples, choose one at random...
          ind_hmm = randi(Kb);
          A = h3m_b.hmm{ind_hmm}.A;
          % perturb A
          A = A + (.5/N) * A;
          % normalize A
          A     = A ./ repmat(sum(A,2),1,N);
          emit = cell(1,N);
          for n = 1 : N
            emit{n} = reduced_out;
            % and randomize priors ...
            emit{n}.priors = rand(size(emit{n}.priors));
            emit{n}.priors = emit{n}.priors / sum(emit{n}.priors);
          end
          
          h3m.hmm{j}.prior = prior;
          h3m.hmm{j}.A     = A;
          h3m.hmm{j}.emit  = emit;
        end
        omega = rand(1,Kr);
        omega = omega/sum(omega);
        h3m.omega = omega;
        h3m.LogL = -inf;
        
      case {'gmm_A'}
        
        h3m.K = Kr;
        ind_hmm = randperm(Kb);
        
        for j = 1 : Kr
          prior = h3m_b.hmm{ind_hmm(j)}.prior;
          % pick one A form the examples, choose one at random...
          
          A = h3m_b.hmm{ind_hmm(j)}.A;
          emit = cell(1,N);
          for n = 1 : N
            emit{n} = reduced_out;
            % and randomize priors ...
            emit{n}.priors = rand(size(emit{n}.priors));
            emit{n}.priors = emit{n}.priors / sum(emit{n}.priors);
          end
          
          h3m.hmm{j}.prior = prior;
          h3m.hmm{j}.A     = A;
          h3m.hmm{j}.emit  = emit;
        end
        omega = rand(1,Kr);
        omega = omega/sum(omega);
        h3m.omega = omega;
        h3m.LogL = -inf;
        
        
      case {'gmm_L2R'}
        
        h3m.K = Kr;
        for j = 1 : Kr
          prior = rand(N,1);
          prior = prior/sum(prior);
          
          % pick one A form the examples, choose one at random...
          ind_hmm = randi(Kb);
          A = h3m_b.hmm{ind_hmm}.A;
          % perturb A
          %A = A + (.5/N) * A;
          % normalize A
          A     = A ./ repmat(sum(A,2),1,N);
          emit = cell(1,N);
          for n = 1 : N
            emit{n} = reduced_out;
            % and randomize priors ...
            emit{n}.priors = rand(size(emit{n}.priors));
            emit{n}.priors = emit{n}.priors / sum(emit{n}.priors);
          end
          
          h3m.hmm{j}.prior = prior;
          h3m.hmm{j}.A     = A;
          h3m.hmm{j}.emit  = emit;
        end
        omega = rand(1,Kr);
        omega = omega/sum(omega);
        h3m.omega = omega;
        h3m.LogL = -inf;
        
        % how many gmms? as many as the number of
        
      case {'g3k'} % a gmm for each hmm component ...
        
        [reduced_out] = cluster_gmms_init(gmms,Kr, virtualSamples, reduced_out,iterations, trials);
        
        h3m.K = Kr;
        for j = 1 : Kr
          prior = rand(N,1);
          prior = prior/sum(prior);
          A     = rand(N);
          A     = A ./ repmat(sum(A,2),1,N);
          emit = cell(1,N);
          for n = 1 : N
            emit{n} = reduced_out.mix{j};
            % and randomize priors ...
            emit{n}.priors = rand(size(emit{n}.priors));
            emit{n}.priors = emit{n}.priors / sum(emit{n}.priors);
          end
          
          h3m.hmm{j}.prior = prior;
          h3m.hmm{j}.A     = A;
          h3m.hmm{j}.emit  = emit;
        end
        omega = rand(1,Kr);
        omega = omega/sum(omega);
        h3m.omega = omega;
        h3m.LogL = -inf;
      case {'g3k_Ad'} % a gmm for each hmm component ...
        
        [reduced_out] = cluster_gmms_init(gmms,Kr, virtualSamples, reduced_out,iterations, trials);
        
        h3m.K = Kr;
        for j = 1 : Kr
          prior = rand(N,1);
          prior = prior/sum(prior);
          
          % pick one A form the examples, choose one at random...
          ind_hmm = randi(Kb);
          A = h3m_b.hmm{ind_hmm}.A;
          % perturb A
          A = A + (.5/N) * A;
          % normalize A
          A     = A ./ repmat(sum(A,2),1,N);
          emit = cell(1,N);
          for n = 1 : N
            emit{n} = reduced_out.mix{j};
            % and randomize priors ...
            emit{n}.priors = rand(size(emit{n}.priors));
            emit{n}.priors = emit{n}.priors / sum(emit{n}.priors);
          end
          
          h3m.hmm{j}.prior = prior;
          h3m.hmm{j}.A     = A;
          h3m.hmm{j}.emit  = emit;
        end
        omega = rand(1,Kr);
        omega = omega/sum(omega);
        h3m.omega = omega;
        h3m.LogL = -inf;
        
        
        
        
      case {'g3n'} % a gmm for each state ...
        
        
        [reduced_out] = cluster_gmms_init(gmms,N, virtualSamples,reduced_out, iterations, trials);
        
        h3m.K = Kr;
        for j = 1 : Kr
          prior = rand(N,1);
          prior = prior/sum(prior);
          A     = rand(N);
          A     = A ./ repmat(sum(A,2),1,N);
          emit = cell(1,N);
          for n = 1 : N
            emit{n} = reduced_out.mix{n};
            % and randomize priors ...
            emit{n}.priors = rand(size(emit{n}.priors));
            emit{n}.priors = emit{n}.priors / sum(emit{n}.priors);
          end
          
          h3m.hmm{j}.prior = prior;
          h3m.hmm{j}.A     = A;
          h3m.hmm{j}.emit  = emit;
        end
        omega = rand(1,Kr);
        omega = omega/sum(omega);
        h3m.omega = omega;
        h3m.LogL = -inf;
        
      case {'g3n_Ad'} % a gmm for each state ...
        
        
        [reduced_out] = cluster_gmms_init(gmms,N, virtualSamples,reduced_out, iterations, trials);
        
        h3m.K = Kr;
        for j = 1 : Kr
          prior = rand(N,1);
          prior = prior/sum(prior);
          
          % pick one A form the examples, choose one at random...
          ind_hmm = randi(Kb);
          A = h3m_b.hmm{ind_hmm}.A;
          % perturb A
          A = A + (.5/N) * A;
          % normalize A
          A     = A ./ repmat(sum(A,2),1,N);
          emit = cell(1,N);
          for n = 1 : N
            emit{n} = reduced_out.mix{n};
            % and randomize priors ...
            emit{n}.priors = rand(size(emit{n}.priors));
            emit{n}.priors = emit{n}.priors / sum(emit{n}.priors);
          end
          
          h3m.hmm{j}.prior = prior;
          h3m.hmm{j}.A     = A;
          h3m.hmm{j}.emit  = emit;
        end
        omega = rand(1,Kr);
        omega = omega/sum(omega);
        h3m.omega = omega;
        h3m.LogL = -inf;
    end
    
                         
end

% debugging - show the initial h3m
if 0
  tmp_h3m = h3m;
  tmp_h3m.Z = [];
  tmp_h3m.LogLs = [];
  
  tmp_hmms = h3m_to_hmms(tmp_h3m);
  figure
  vhem_plot(tmp_hmms, 'face.jpg');
  drawnow
end

% make the prior and A
function [prior, A] = makeAprior(N, mode)
if any(mode == 'u')
  prior = ones(N,1)/N;
  A     = ones(N,N)/N;
elseif any(mode == 'r')
  prior = rand(N,1);
  prior = prior/sum(prior);
  A = rand(N);
  A = A ./ repmat(sum(A,2),1,N);
else
  error('unknown mode');
end

% make the GMM component weights
function [p_out] = makeGMMweights(p, A, mode)
p = p(:)';
  
if (any(mode=='0'))
  for t = 1 : 50
    p = p * A;
  end
  p_out = p;
  
elseif (any(mode=='3'))
  T = 3;
  p_out = p;
  for t=2:T
    p = p * A;
    p_out = p_out + p;
  end  
  p_out = p_out/T;

elseif (any(mode=='x'))
  p_out = ones(size(p))/size(p,2);

else
  error('unknown mode');
end

%% plot a Gaussian as ellipse
function plot2D(mu, Sigma, color)
% truncate to 2D
mu = mu(1:2);
Sigma = Sigma(1:2,1:2);

mu = mu(:);
if ~any(isnan(Sigma(:))) && ~any(isinf(Sigma(:)))
  [U,D] = eig(Sigma);
  n = 100;
  t = linspace(0,2*pi,n);
  xy = [cos(t);sin(t)];
  k = sqrt(conf2mahal(0.95,2));
  w = (k*U*sqrt(D))*xy;
  z = repmat(mu,[1 n])+w;
  h = plot(z(1,:),z(2,:),color,'LineWidth',1);
end

function m = conf2mahal(c,d)
m = chi2inv(c,d);