function emit = hem_fix_degenerate_emission(emit, i_z)
% hem_fix_degenerate_emission - fix a degenerate component of a GMM emission [internal function]
%
% emit = hem_fix_degenerate_emission(emit, i_z)
% 
% i_z is the index to of the GMM emission to fix
%
% ---
% H3M Toolbox 
% Copyright 2012, Emanuele Coviello, CALab, UCSD
% Antoni Chan 2017/06/06

fprintf('!!! modifying gmm emission: one component is zero \n')

[foo highest] = max(emit.priors);
emit.priors([i_z highest]) = emit.priors(highest)/2;

% renormalize for safety
emit.priors = emit.priors / sum(emit.priors);
                
emit.centres(i_z,:) = emit.centres(highest,:);
switch(emit.covar_type)
  case 'diag'
    emit.covars(i_z,:)  = emit.covars(highest,:);
  case 'full'
    error('here');
end

                
% perturb only the centres
emit.centres(i_z,:) =  emit.centres(i_z,:) + (0.01 * rand(size(emit.centres(i_z,:))))...
                     .* emit.centres(i_z,:);
