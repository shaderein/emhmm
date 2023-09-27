function str = prob2str(p, numplaces)
% prob2str - convert probability to a string
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-01-13
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong

if nargin<2
  numplaces = 2;
end

str = sprintf(sprintf('%%.%df', numplaces), p);
if str(1) == '0'
  str = str(2:end);
else
  % it should be 1.0
  str = str(1:end-1);
end

