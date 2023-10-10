function vbco = vbhem_co_permute_clusters(vbco, cl)
%  vbhem_co_permute_clusters - permute the groups in a VBHEM co-clustering
% 
%    group_hmms_new = vbhem_co_permute_clusters(group_hmms, cl)
% 
%  maps group "cl(i)" in group_hmms to group "i" in group_hmms_nmew

% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2021-09-02
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

% 2021-07-29: v0.80 - ABC - initial version

% swap HMMs 
for i=1:length(vbco.cogroup_hmms)
  tmp = vbco.cogroup_hmms{i};
  
  % swap hmms
  tmp_new = tmp;
  tmp_new.alpha      = tmp.alpha(cl);
  tmp_new.omega      = tmp.omega(cl);
  tmp_new.Z          = tmp.Z(:,cl);
  tmp_new.Nj         = tmp.Nj(cl);
  %tmp_new.L_elbo     = tmp.L_elbo(:,:,cl);
  tmp_new.label      = cl(tmp.label);
  tmp_new.groups     = tmp.groups(cl);
  tmp_new.group_size = tmp.group_size(cl);
  tmp_new.hmms       = tmp.hmms(cl);
    
  vbco.cogroup_hmms{i} = tmp_new;

  % error check
  handledfields = {'alpha', 'omega', 'Z', 'Nj', 'label', 'groups', ...
    'group_size', 'hmms', 'K', 'S', 'W0mode'};
  error_check(handledfields, tmp_new);
end

% copy the group info
vbco.label      = vbco.cogroup_hmms{1}.label;
vbco.groups     = vbco.cogroup_hmms{1}.groups;
vbco.group_size = vbco.cogroup_hmms{1}.group_size;
vbco.NL         = vbco.cogroup_hmms{1}.Nj;

% swap hmms
vbco.L_elbo     = vbco.L_elbo(:,:,cl);

% error check
handledfields = {'K', 'label', 'groups', 'group_size', 'iter', 'LL', 'L_elbo', ...
  'NL', 'learn_hyps', 'model_LL', 'model_k', 'model_bestK', 'model_all', ...
  'cogroup_hmms', 'vbhemopt', 'emhmm_version', 'matlab_version'};
error_check(handledfields, vbco);


% check unhandled fields %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function error_check(handledfields, A)

tmp = fields(A);
for i=1:length(handledfields)
  foo = ~strcmp(handledfields{i}, tmp);
  tmp = tmp(foo);
end

% issue warning if some fields are not handled
if ~isempty(tmp)
  tmp
  error('need to update this function!');
end




