function [h3mo, emptyinds] = vbh3m_remove_empty(h3m, thresh,verbose)
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2021-07-29
% Antoni B. Chan, Janet H. Hsiao, Lan Hui
% City University of Hong Kong, University of Hong Kong

% 2021-07-29: v0.80 - initial version

if (nargin<2)
  thresh = 1.0;
end

if (nargin<3)
  verbose = 1;
end
%if (nargin<4)
%  sortclusters= 'f';
%end

N = h3m.Nj;
check  = (N < thresh);
zinds  = find(check);
nzinds = find(~check);

h3mo = h3m;

% remove empty hmms and update information
if ~isempty(zinds)
  if (verbose)
    fprintf('Removing hmms: %s --> K=%d\n', mat2str(zinds(:)'), length(nzinds));
  end
  
  
  h3mo.hmms = h3m.hmms(nzinds);
  h3mo.omega = h3m.omega(nzinds);
  % 2020-9-10: nomorlize it 
  h3mo.omega = h3mo.omega./sum(h3mo.omega);
  h3mo.Z = h3m.Z(:,nzinds);
  if isfield(h3mo, 'NL')
    h3mo.NL = h3m.NL(:,nzinds);
  end
  if isfield(h3mo, 'L_elbo')
    h3mo.L_elbo = h3m.L_elbo(:,nzinds);
  end
  for i = 1:length(nzinds)
    h3mo.label(h3m.label== nzinds(i))=i;
  end
  h3mo.K = length(h3mo.hmms);
  

  % update variation parameters

%   h3mo.varpar.eta   = h3m.varpar.eta(nzinds);
%   h3mo.varpar.alpha   = h3m.varpar.alpha(nzinds);
%   h3mo.varpar.epsilon = h3m.varpar.epsilon(nzinds);
% 
%   h3mo.varpar.beta    = h3m.varpar.beta(:,:,nzinds);
%   h3mo.varpar.v       = h3m.varpar.v(:,:,nzinds);
%   h3mo.varpar.m       = h3m.varpar.m(nzinds);
%   h3mo.varpar.W       = h3m.varpar.W(nzinds);
  
  % update counts
%   h3mo.NJ_rho1 = h3m.NJ_rho1 (:,nzinds);
%   h3mo.NJ_ror = h3m.NJ_ror(:,:, nzinds);
%   h3mo.NJ_rho = h3m.NJ_rho(nzinds,:);
%   h3mo.NJ  = h3m.NJ(1,nzinds);
  
  % update groups
  h3mo.groups = h3m.groups(1,nzinds);
  h3mo.group_size = h3m.group_size(1,nzinds);

  emptyinds = zinds;
else
  % do nothing
  emptyinds = [];
  
end

% remove empty states
hmms= h3mo.hmms;

for i=1:length(hmms)
    [hmms{i}, zi] = vbhmm_remove_empty(hmms{i}, 0, 1e-3);
    if ~isempty(zi)
      if (verbose)
        fprintf('HMM %d: removed states', i);
        fprintf(' %d', zi);
        fprintf(' --> S=%d\n', length(hmms{i}.pdf));
      end
    end

end
h3mo.hmms = hmms;

% sort the ROIs (since some may be removed) 
%for j=1:length(h3mo.hmms)
%    hmms = h3mo.hmms{j};
%    h3mo.hmms{j} = vbhmm_standardize(hmms, sortclusters);
%end
