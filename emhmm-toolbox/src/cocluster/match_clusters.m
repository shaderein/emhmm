function [newinds] = match_clusters(L_elbo, mode, opt, showfig)
% match_clusters - match clusters betweeen HMMs for initializing co-clustering (internal function)
%
% L_elbo = [stimuli x Kb x Kr]
% mode = 'r' -- random select starting index, and randomly select best matches
% mode = 'g' -- greedy select best matches, starting at index "opt"
% mode = 's' -- spectral clustering, opt is an annealing parameter
% mode = 'k' -- k-means clustering
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2019-07-13
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

% 2020-03-12: v0.77 - initial version


%
%          h3m.hmm = indh3m{i}.hmm(newinds(i,:));
%
% newinds(stimuli, :) = [index permutation list]

if nargin<2
  mode = 0;
end
if nargin<3
  showfig = 0;
end

[NS, Kb, Kr] = size(L_elbo);

% compute posterior for each stimuli
% (ignore virtual samples and omega)
log_Z = zeros(NS,Kb,Kr);
Z     = zeros(NS,Kb,Kr);
for k=1:NS
  tmp = reshape(L_elbo(k,:,:),[Kb,Kr]);  
  log_Z(k,:,:) = tmp - logtrick(tmp')'*ones(1,Kr);
end
Z = exp(log_Z);

% [stimuli x Kb x Kr]
%SCORE = L_elbo;
SCORE = Z;


% permutations for Kr
allperms = perms(1:Kr);
numperms = size(allperms,1);

switch(mode)
  
  %% RANDOM order
  case 'r'
    
    % random order
    order = randperm(NS);
    myk = order(1);
    
    % use the first stimuli as the base
    newinds = zeros(NS,Kr);
    newinds(myk,:) = 1:Kr;
    
    % current score [Kb x Kr]
    score = squeeze(SCORE(myk,:,:));
    
    % go through each other stimuli
    for k=2:NS
      if showfig
        figure(111)
        clf
      end
      
      myk = order(k);
      
      bestmax = -inf;
      bestperm = [];
      bestscore = [];
      bestind = [];
      for j=1:numperms
        % compute new score using this permutation
        myperm = allperms(j,:);
        newscore = score + squeeze(SCORE(myk, :, myperm));
        
        % get max over each row and sum to get consistency score
        newmax = sum(max(newscore, [], 2));
        
        % save if it's better
        if newmax>bestmax
          bestmax = newmax;
          bestperm = myperm;
          bestscore = newscore;
          bestind = j;
        end
        
        if showfig
          subplot(1,numperms+1,j+1);
          imagesc(newscore);
          title(sprintf('%g', newmax));
        end
      end
      
      if showfig
        subplot(1,numperms+1,1)
        imagesc(score);
        title('original');
        subplot(1,numperms+1,bestind+1)
        title(sprintf('%g best', bestmax));
        colormap(gray)
        drawnow
        %pause
        
      end
      
      % save best permutation
      score = bestscore;
      newinds(myk,:) = bestperm;  
    end
    
    
  %% fixed greedy order  
  case 'g'
    myk = opt;
    
    % use the first stimuli as the base
    newinds = zeros(NS,Kr);
    newinds(myk,:) = 1:Kr;
    
    % current score [Kb x Kr]
    score = squeeze(SCORE(myk,:,:));
    
    usedk = zeros(1,NS);
    usedk(myk) = 1;
    
    % go through each other stimuli
    for k=2:NS
      if showfig
        figure(111)
        clf
      end
      
      bestmax = -inf;
      bestperm = [];
      bestscore = [];
      bestind = [];
      bestk   = [];
      
      % find best permutation among all remaining stimuli
      for myk=1:NS
        if ~usedk(myk)
          for j=1:numperms
            % compute new score using this permutation
            myperm = allperms(j,:);
            newscore = score + squeeze(SCORE(myk, :, myperm));
            
            % get max over each row and sum to get consistency score
            newmax = sum(max(newscore, [], 2));
            
            % save if it's better
            if newmax>bestmax
              bestmax = newmax;
              bestperm = myperm;
              bestscore = newscore;
              bestind = j;
              bestk = myk;
            end
          end
        end
      end
      
      if showfig
        subplot(1,2,1)
        imagesc(score);
        title('original');
        subplot(1,2,2)
        imagesc(bestscore)
        colormap(gray)
        title(sprintf('best %d: %g', bestk, bestmax));
        drawnow
        %pause
      end
      
      % save best permutation
      score = bestscore;
      newinds(bestk,:) = bestperm;
      usedk(bestk) = 1;
    end

  case 's'
    %% Spectral Clustering
    
    % initialize probability matrix
    P = NS*eye(Kb, Kb);
    
    % for each stimuli
    for k=1:NS
      % for each pair (q,r)
      for q=1:Kb
        for r=(q+1):Kb          
          % compute probability that q and r are assigned to the same class
          qcl = squeeze(SCORE(k,q,:));
          rcl = squeeze(SCORE(k,r,:));
          tmp = sum(qcl.*rcl);
          P(q,r) = P(q,r)+tmp;
          P(r,q) = P(q,r);
        end
      end      
    end
    
    % normalize
    P = P/NS;
    
    % annealing parameter
    P = P.^opt;
    
    
    % run spectral clustering 
    [ C ] = SpectralClustering(P, Kr, Kr, 1);    
    C = full(C);
    
    % get class assignments
    [~,cl] = max(C, [], 2);       
    
    newinds = zeros(NS,Kr);
        
    % for each stimuli, pick the permutation that is most similar to C
    for k=1:NS      
      bestscore = -inf;
      bestinds  = [];
      for j=1:numperms
        % compute new score using this permutation
        myperm = allperms(j,:);
        
        tmp = squeeze(SCORE(k,:,myperm));
        tmpscore = sum(sum(tmp.*C));
        if (tmpscore>bestscore)
          bestinds = myperm;
          bestscore = tmpscore;
        end        
      end
      
      newinds(k,:) = bestinds;      
    end
    
    
    if (showfig)
      % reorder index
      inds = [];
      for k=1:Kr
        inds = [inds; find(cl==k)];
      end
      
      subplot(1,3,1)
      imagesc(P)
      subplot(1,3,2)
      imagesc(P(inds,inds));
      subplot(1,3,3)
      
      allscore = zeros(NS+2, Kb*Kr);
      allscore(1,:) = C(:);
      for k=1:NS
        tmp = SCORE(k,:,newinds(k,:));
        allscore(k+1,:) = tmp(:);
      end
      
      allscore(end,:) = mean(allscore(2:end-1,:), 1);
      imagesc(allscore);
      
    end
    
  case 'k'
    %% use k-means clustering on score vectors
  
    %SCORE
 
    % [NS x Kb x Kr] -> [Kb x NS x Kr]
    X = permute(SCORE, [2 1 3]);
    
    X2 = reshape(X, [size(X,1), size(X,2)*size(X,3)]);
        
    % run K-means
    [idx, c] = kmeans(X2, Kr, 'Start', 'sample');
    
    % [Kr x (NS*Kr)] -> [Kr x NS x Kr]
    cs = reshape(c, [Kr, size(X,2), size(X,3)]);
    
    % for each stimuli
    for i=1:NS
      A = squeeze(cs(:,i,:));
      myscores = zeros(1,numperms);
    
      for n=1:numperms
        myperm = allperms(n,:);
        % rearrange and then sum the diagonal
        Atmp = A(:, myperm);
        myscores(n) = trace(Atmp);
      end
      
      [~,maxind] = max(myscores);
      
      newinds(i,:) = allperms(maxind, :);
    end
    
    
    
    if (showfig)
      tmp = zeros(1,NS*Kr);
      for i=1:NS
        tmp(i+(0:(Kr-1))*NS) = (newinds(i,:)-1)*NS + i;
      end
      %newinds
      %tmp
        
      
      subplot(3,2,1)
      imagesc(X2);
      for j=1:2
        ii = find(idx==j);
        subplot(3,2,1+2*j)
        imagesc([c(j,:); X2(ii,:)])
        
        subplot(3,2,2+2*j);        
        imagesc([c(j,tmp); X2(ii,tmp)]);        
      end
    
      
      
    end
    
end


