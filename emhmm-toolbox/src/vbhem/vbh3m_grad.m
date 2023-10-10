function [L, dL] = vbh3m_grad(h3m_b, vbopt, hypinfo)
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2021-07-29
% Antoni B. Chan, Janet H. Hsiao, Lan Hui
% City University of Hong Kong, University of Hong Kong

% 2021-07-29: v0.80 - initial version


% calculate derivative?
if (nargout>1)
  vbopt.calc_LLderiv = 1;
end

% % run EM
[vbopt, hyp_clipped] = vbhem_clip_hyps(vbopt);
vbopt.hyp_clipped = hyp_clipped;
%             
% initialization
[h3m_r,vbopt] = vbhemhmm_init(h3m_b,vbopt);

h3m = vbhem_h3m_c_step_fc(h3m_r,h3m_b,vbopt);

L = - h3m.LL;

if (nargout>1)  
  dL = [];
  % append gradients
  for i=1:length(hypinfo)
    dL = [dL; -h3m.dLL.(hypinfo(i).derivname)(:)];
  end
end


if (vbopt.verbose >=2)
  fprintf('** gradient: dL=');
  fprintf('(');
  fprintf('%g ', -dL);
  fprintf(')\n');
end
