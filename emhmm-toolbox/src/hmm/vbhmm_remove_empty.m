function [hmmo, emptyinds] = vbhmm_remove_empty(hmm, verbose, thresh)
% vbhmm_remove_empty - remove empty states in an HMM
%
%   [hmmo, emptyinds] = vbhmm_remove_empty(hmm, verbose, thresh)
%
%    hmm = input HMM
%    verbose = 0 : quiet
%              1 : verbose [default]
%    thresh = threshold on N for removing empty states 
%             [default = 1.0, at least 1 fixation per ROI]
%
%    hmmo = output HMM - empty states are removed
%   emptyinds = states that were removed, [] if none
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2018-02-13
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

% 2018-03    - v0.72 initial version
% 2018-04-20 - v0.72 support for groups
% 2018-05-15 - v0.72 settable threshold, and set to 1.0. fix a few bugs (N, normalize gamma)
% 2018-06-15 - v0.74 add emptyinds output       
% 2021-07-29 - v0.80 handle no gamma (Hui)

if (nargin<2)
  verbose = 1;
end
if (nargin<3)
  thresh = 1.0;
end

N = hmm.N;
check  = (N < thresh);
zinds  = find(check);
nzinds = find(~check);

hmmo = hmm;

usegroups = isfield(hmm, 'group_map');
if usegroups
  G = length(hmm.group_ids);
end

if ~isempty(zinds)
  if (verbose)
    fprintf('removing states: %s\n', mat2str(zinds(:)'));
  end
  
  % update variation parameters
  if ~usegroups
    hmmo.varpar.alpha   = hmm.varpar.alpha(nzinds);
    hmmo.varpar.epsilon = hmm.varpar.epsilon(nzinds, nzinds);
  else
    for g=1:G
      hmmo.varpar.alpha{g}   = hmm.varpar.alpha{g}(nzinds);
      hmmo.varpar.epsilon{g} = hmm.varpar.epsilon{g}(nzinds, nzinds);
    end
  end
  hmmo.varpar.beta    = hmm.varpar.beta(nzinds);
  hmmo.varpar.v       = hmm.varpar.v(nzinds);
  hmmo.varpar.m       = hmm.varpar.m(:,nzinds);
  hmmo.varpar.W       = hmm.varpar.W(:,:,nzinds);
  
  % update counts
  if ~usegroups
    hmmo.M  = hmm.M(nzinds, nzinds);
    hmmo.N1 = hmm.N1(nzinds);
  else
    for g=1:G
      hmmo.M{g}  = hmm.M{g}(nzinds, nzinds);
      hmmo.N1{g} = hmm.N1{g}(nzinds);
      hmmo.Ng{g} = hmm.Ng{g}(nzinds);
    end
  end
  hmmo.N  = hmm.N(nzinds);
  
  if isfield(hmm, 'gamma')
    for i=1:length(hmm.gamma)
      hmmo.gamma{i} = hmm.gamma{i}(nzinds,:);
      % renormalize gamma
      hmmo.gamma{i} = bsxfun(@times, hmmo.gamma{i}, 1./sum(hmmo.gamma{i}, 1));
    end
  end
  
  % update parameters  
  if ~usegroups
    hmmo.prior = hmmo.varpar.alpha ./ sum(hmmo.varpar.alpha);
    hmmo.trans = hmmo.varpar.epsilon ./ ...
      repmat(sum(hmmo.varpar.epsilon, 2), [1,length(nzinds)]);
  else
    for g=1:G
      hmmo.prior{g} = hmmo.varpar.alpha{g} ./ sum(hmmo.varpar.alpha{g});
      hmmo.trans{g} = hmmo.varpar.epsilon{g} ./ ...
        repmat(sum(hmmo.varpar.epsilon{g}, 2), [1,length(nzinds)]);
    end
  end
  hmmo.pdf   = hmm.pdf(nzinds);
  
  emptyinds = zinds;
else
  % do nothing
  emptyinds = [];
  
end