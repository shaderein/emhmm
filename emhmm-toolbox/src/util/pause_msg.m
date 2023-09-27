function pause_msg(do_pause, msg)
% pause_msg - pause and show a message
%
%  pause_msg(do_pause, msg)
%
%   do_pause = 1 pause, 0 don't pause
%        msg = message to show if pausing
%              default = '<press a key to continue>'
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2019-06-14
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

% v0.76 - initial version

if nargin<2
  msg = '<press any key to continue>';
end

if (do_pause)
  fprintf(['<strong>' msg '</strong>\n']);
  pause
end
