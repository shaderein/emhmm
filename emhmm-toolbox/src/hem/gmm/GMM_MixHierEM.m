function [mix, post, logpost, redo] = GMM_MixHierEM(allMix,T, virtualSamples, iterations)
%
%Implemenation of Mixture Hierarchies Dessity Estimation
% EM for combining mixture of Gaussian Distrubtions into one distribution
% [Vasconcelos '01, Carneiro Vasconcelos '05]. 
% NOTE: Upper level centers that are not responsible for any datapoints are
% randomly reassigned.
% NOTE: Lower level centers that have zero likelihood under all layer
% components are ignored.
%
%
% Implemented by Douglas Turnbull  - March 2006
%
% INPUTS:
%  allMix is an array of gaussian mixture all with K components
%  T number of components in resulting mixture (default - size of allMix)
%  virtualSamples - positive scalar number of samples to be drawn from
%    each allMix mixtures (default - 1 per GMM in all Mix).
%  createMovie = boolean 1 to create a movie of mixture distribution, 0
%  otherwise - M will be []
%
% OUTPUTS:
%  mix is the resulting mixture after using mixture Hiearchies EM
%
% Need to implement 
%  1. Sanjoy's extra centers/pruning algorithm
%  2. Sanjoy's initialization using furthest first variant

% 2017-08-16: ABC: added full covariance support
% 2018-06-08: v0.73 - tested full covariance support for T=1
% 2018-06-15: v0.74 - improve numerical stability for posterior

redo = 0;
if nargin < 2, T = max(size(allMix)); end
if nargin < 3, virtualSamples = length(allMix)*allMix(1).mix.ncentres; end
if nargin < 4, iterations = 100; end

diffThreshold = 0.0010; %theshold at which if average percentage change in
                        % mix.priors, mix.centres, and mix.covars for all
                        % parameters is less, then we stop parameter
                        % estimation.

%number of song distributions
D = max(size(allMix));
dim= allMix(1).mix.nin;

%number of components for each song distribution must be the same.
K = zeros(D,1);
for j = 1:D, K(j) = allMix(j).mix.ncentres; end
if min(K) ~= max(K)
    error('hierarchyWordMix: All mixture distributions must have same number of components.');
    return
end
K = max(K);

covar_type = allMix(1).mix.covar_type;

%Combine LOWER LEVEL MIXTURE
naiveMix = naiveWordMix(allMix);

data = naiveMix.centres';  % [dim x npts]
switch(covar_type)
  case 'diag'
    dataCov = naiveMix.covars';  % [dim x npts]
  case 'full'
    dataCov = naiveMix.covars;   % [dim x dim x npts]
end
dataPrior = repmat(naiveMix.priors,T,1);
[ndim,npts] = size(data);
%
if T ==1
  
  
  mix = gmm(ndim,T, covar_type); %initialize structure
  mix.priors  = 1;
  mix.centres = dataPrior*data';
  switch(covar_type)
    case 'diag'
      mix.covars  = dataPrior*((data').^2+dataCov');
    case 'full'
      data2 = bsxfun(@times, reshape(data,[dim 1 npts]), reshape(data,[1 dim npts]));
      mix.covars = sum(bsxfun(@times, reshape(naiveMix.priors, [1 1 npts]), ...
        data2 + dataCov), 3);
  end
  
  post = ones(T,size(data,2));
  logpost = zeros(size(post));  
  
else
  
  post = rand(T,size(data,2));
  
  cent = arKmeans(data', T,1,0)';
  switch(covar_type)
    case 'diag'  % [dim x T]
      vrnc = repmat(max(naiveMix.covars), T,1)'; % ALSO TRY MEAN instead of MAX
    case 'full'  % [dim x dim x T]
      vrnc = repmat(mean(naiveMix.covars,3), [1 1 T]); 
  end
  
  mxwt = ones(T,1)/T;
  %preallocate arrays to check for early stopping
  % measures average difference in the value of each parameter between iterations
  coef = -(ndim/2)*log(2*pi);
  dpp = dataPrior .*virtualSamples;
  
  %EM Stopping parameter
  last = -realmax;
  tol = 1e-6;
  
  for  i = 1:iterations
    
    %E-step
    switch(covar_type)
      case 'diag'
        ivr = 1./vrnc;
        trCov = ivr'*dataCov; % [T npts]
    
        nrm = sum(cent.*cent.*ivr + log(vrnc),1)';
        %xpt is the log of the numerator of the posterior. [T npts]
        xpt = log(mxwt(:,ones(1,npts))) + dpp .* (coef ...
          - 0.5*(ivr'*(data.*data) - 2*(cent.*ivr)'*data + nrm(:,ones(1,npts))+trCov));

      case 'full'
        ivr = zeros(dim,dim,T);
        ld_vrnc = zeros(T,1);
        trCov = zeros(T,npts);
        for t=1:T
          ivr(:,:,t) = inv(vrnc(:,:,t));
          tmp = bsxfun(@times, ivr(:,:,t), dataCov); 
          trCov(t,:) = sum(sum(tmp,1), 2);
          
          ld_vrnc(t) = log(det(vrnc(:,:,t)));
        end
        
        tmp = zeros(T, npts);
        for t=1:T
          for n=1:npts
            foo = cent(:,t) - data(:,n);
            tmp(t,n) = foo'*ivr(:,:,t)*foo;
          end
        end
                
        xpt = log(mxwt(:,ones(1,npts))) + ...
          dpp .* (coef - 0.5*(trCov+tmp+repmat(ld_vrnc,[1 npts])));
        
    end
    
    % POSTERIOR
    if 0
      maxx = max(xpt);
      post = exp(xpt-repmat(maxx,T,1));
      prob = sum(post,1);
      post = post./repmat(prob,T,1);
      logp = mean(log(prob+realmin)+maxx); %/size(data,1);
    end
    
    % 2018-06-15 - v0.74 - use more stable computations
    lxpt = logtrick(xpt);
    % ABC - log posterior   
    logpost = bsxfun(@minus, xpt, lxpt);
    % posterior
    post = exp(logpost);
    % logp
    logp = mean(lxpt);
    
    
    if (~isfinite(logp))
      redo = 1;
      fprintf('re-doing the HEM-GMM')
      break
      %         error('Infinity in EM');
    end;
    
    %disp(sprintf('Mixture Hierarchies EM: %d\tlogp: %f',i,logp));
    if ((logp-last)<tol), break; end;
    last = logp;
    
    %M-step
    mxwt = mean(post',1)';
    wts =  post.*dataPrior;
    wts =  wts./ repmat(sum(wts,2),1,npts);
    cent = data*wts';
    
    %REMOVE THIS LOOP BY CROSS MULTIPLYING - might now work because of
    %outer product
    for c = 1:T
      diffs = data - repmat(cent(:,c),1,npts);
      switch(covar_type)
        case 'diag'
          vrnc(:,c) = (diffs.^2+dataCov)*wts(c,:)';
        case 'full'
          tmp = bsxfun(@times, reshape(diffs, [dim 1 npts]), reshape(diffs, [1 dim npts]));
          vrnc(:,:,c) = sum(bsxfun(@times, tmp+dataCov, reshape(wts(c,:), [1 1 npts])), 3);
          
          %vrnc(:,:,c) = diag(diag(vrnc(:,:,c))); % test diagonal full
      end
    end
    
  end
  
  
  mix = gmm(ndim,T, covar_type); %initialize structure
  mix.priors  = mxwt';
  mix.centres = cent';
  switch(covar_type)
    case 'diag'
      mix.covars  = vrnc';
    case 'full'
      mix.covars = vrnc;
  end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s] = logtrick(lA)
% logtrick - "log sum trick" - calculate log(sum(A)) using only log(A) 
%
%   s = logtrick(lA)
%
%   lA = column vector of log values
%
%   if lA is a matrix, then the log sum is calculated over each column
% 

[mv, mi] = max(lA, [], 1);
temp = lA - repmat(mv, size(lA,1), 1);
cterm = sum(exp(temp),1);
s = mv + log(cterm);




