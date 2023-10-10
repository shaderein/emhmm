function group_hmms_new = vhem_permute_clusters(group_hmms, cl)
% vhem_permute_clusters - permute the groups in a VHEM clustering
%
%   group_hmms_new = vhem_permute_clusters(group_hmms, cl)
%
% maps group "cl(i)" in group_hmms to group "i" in group_hmms_nmew
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2021-07-20
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong

% v0.78 - initial version

%
% maps group for cl(i)-> i
%

group_hmms_new = group_hmms;

group_hmms_new.Z          = group_hmms.Z(:,cl);
group_hmms_new.label      = cl(group_hmms.label);
group_hmms_new.groups     = group_hmms.groups(cl);
group_hmms_new.group_size = group_hmms.group_size(cl);
group_hmms_new.hmms       = group_hmms.hmms(cl);


handledfields = {'Z', 'label', 'groups', 'group_size', 'hmms', ...
  'LogLs', 'LogL', 'hemopt', 'emhmm_version', 'matlab_version'};
tmp = fields(group_hmms);
for i=1:length(handledfields)
  foo = ~strcmp(handledfields{i}, tmp);
  tmp = tmp(foo);
end

% issue warning if some fields are not handled
if ~isempty(tmp)
  for i=1:length(tmp)
    warning('field "%s" was not handled. need to update this function.', tmp{i});
  end
end

