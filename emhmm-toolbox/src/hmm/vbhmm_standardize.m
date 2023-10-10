function hmm_new = vbhmm_standardize(hmm, mode)
% vbhmm_standardize - standardize an HMM's states (ROIs) to be consistent.
%
%   hmm_new = vbhmm_standardize(hmm, mode)
%
%     hmm = an HMM from vbhmm_learn
%           or a group HMM from vhem_cluster
%
%    mode = 'e' - sort by emission frequency (overall number of fixations in an ROI)
%         = 'p' - sort by prior frequency (number of first-fixations in an ROI)
%         = 'f' - sort by most-likely fixation path
%                 (state 1 is most likely first fixation. State 2 is most likely 2nd fixation, etc)
%         = 's' - sort by steady-state probability of states
%         = 'l' - sort by emission location (left-to-right)
%         = 'r' - sort by emission location (right-to-left)
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-01-13
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong

% 2016-11-22: ABC - created
% 2018-03-11: v0.72 - support for groups
% 2018-06-08: v0.73 - added sort by steady-state
% 2019-03-08: v0.76 - sort by emission location (left-to-right)

% TODO: sort by emission Gaussian spatial size

% run on each hmm
if isfield(hmm, 'hmms')
  hmm_new = hmm;
  for i=1:length(hmm.hmms)
    hmm_new.hmms{i} = vbhmm_standardize(hmm.hmms{i}, mode);
  end
  
  return
end

switch(mode)
  % reorder based on cluster size
  case {'d', 'e'}
    if mode=='d'
      warning('standardization mode ''d'' is deprecated. Use ''e''');
    end
    
    if isfield(hmm, 'N')
      [wsort, wi] = sort(hmm.N, 1, 'descend');
    else
      error('cluster size unknown - not from vbhmm_learn');
    end
    
    % sort by steady-state probability
  case 's'
    p = vbhmm_prob_steadystate(hmm);
    [wsort, wi] = sort(p, 1, 'descend');
    
  % sort by prior
  case 'p'
    if ~iscell(hmm.prior)
      myprior = hmm.prior;
    else
      % groups: add priors
      myprior = size(hmm.prior{1});
      for g=1:length(hmm.prior)
        myprior = myprior + hmm.prior{g};
      end      
    end    
    [wsort, wi] = sort(myprior(:), 1, 'descend');

  % sort by likely fixation path
  case 'f'
    if ~iscell(hmm.prior)
      A = hmm.trans;
      myp = hmm.prior;
    else
      % groups: use first group
      A = hmm.trans{1};
      myp = hmm.prior{1};
    end
    % find starting point
    for t=1:length(myp)
      % get next most-likely fixation
      if (t==1)
        [~, curf] = max(myp);
      else
        [~, curf] = max(A(curf,:));
      end        
      wi(t) = curf;
      
      % invalidate this fixation
      A(:,curf) = -1;
    end
    
    % sort by emission location (left-to-right, or right-to-left)
  case {'l', 'r'}
    if (mode == 'l')
      m = 'ascend';
    else
      m = 'descend';
    end
    tmp = cellfun(@(x) x.mean(1), hmm.pdf);
    
    [wsort, wi] = sort(tmp(:), 1, m);
    
  otherwise
    error('unknown mode');
end


% permute states
hmm_new = vbhmm_permute(hmm, wi);