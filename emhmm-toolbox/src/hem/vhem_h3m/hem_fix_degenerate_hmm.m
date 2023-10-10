function hmm = hem_fix_degenerate_hmm(hmm, i_z)
% hem_fix_degenerate_hmm - fix a degenerate hmm state [internal function]
%
% hmm = hem_fix_degenerate_hmm(hmm, i_z)
% 
% i_z is the index of the HMM state to fix
%
% ---
% H3M Toolbox 
% Copyright 2012, Emanuele Coviello, CALab, UCSD
% Antoni Chan 2018/02/13

[foo, highest] = max(hmm.stats.emit_vcounts);

hmm.prior([i_z highest]) = hmm.prior(highest)/2;

% renormalize for safety
hmm.prior = hmm.prior / sum(hmm.prior);
                

hmm.A(i_z,:) = hmm.A(highest,:);
hmm.A(:,[i_z highest]) = repmat(hmm.A(:,highest)/2, [1 2]);

hmm.emit{i_z} = hmm.emit{highest};
                
% perturb only the centres
hmm.emit{i_z}.centres =  hmm.emit{i_z}.centres + ...
  (0.01 * rand(size(hmm.emit{i_z}.centres))).* hmm.emit{i_z}.centres;
