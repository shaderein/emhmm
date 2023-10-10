function vbhmm_print_hyps(names, vbopt)
% vbhmm_print_hyps - print hyps
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-08-21
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

fprintf('{');
for k=1:length(names)
  fprintf('%s=%s; ', names{k}, mat2str(vbopt.(names{k}), 4));
end
fprintf('}');

