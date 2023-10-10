function [L, dL] = vbh3m_cograd(h3m_bs, vbopt, hypinfo)
% internal function
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2021-07-29
% Antoni B. Chan, Janet H. Hsiao, Lan Hui
% City University of Hong Kong, University of Hong Kong

% 2021-07-29: v0.80 - Hui - initial version


% calculate derivative?
if (nargout>1)
  vbopt.calc_LLderiv = 1;
end

% run EM
[vbopt, hyp_clipped] = vbhem_clip_hyps(vbopt);
vbopt.hyp_clipped = hyp_clipped;
            
%initialize all stimuli
h3ms =[];
moptpointer = vbopt;
for i=1:vbopt.stimulisize
    if length(vbopt.S)>1
        moptpointer.S = vbopt.S(i);
    end
    moptpointer.inithmm = vbopt.inithmm.h3m_rs{i};
    h3m = vbhemhmm_init(h3m_bs(i),moptpointer);
    h3ms = [h3ms; h3m];
end

h3m_rs = vbhem_coh3m_c_step_fc(h3ms, h3m_bs, vbopt);

%clear h3ms

L = -h3m_rs.LL;

if (nargout>1)  
  dL = [];
  % append gradients
  for i=1:length(hypinfo)
    dL = [dL; -h3m_rs.dLL.(hypinfo(i).derivname)(:)];
  end
end


if (vbopt.verbose >=2)
  fprintf('** gradient: dL=');
  fprintf('(');
  fprintf('%g ', -dL);
  fprintf(')\n');
end
