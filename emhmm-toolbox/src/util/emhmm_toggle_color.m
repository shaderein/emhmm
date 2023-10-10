function emhmm_toggle_color(mode)
% emhmm_toggle_color - toggle between color and grayscale plotting
%
% USAGE:
%   emhmm_toggle_color()  - toggle grayscale/color plotting
%   emhmm_toggle_color(0) - set grayscale plotting
%   emhmm_toggle_color(1) - set color plotting (default)
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2019-02-21
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

% v0.75 - initial version

global EMHMM_COLORPLOTS

names = {'grayscale', 'color'};

if nargin<1
  EMHMM_COLORPLOTS = ~EMHMM_COLORPLOTS;
  
else
  if (mode == 0)
    EMHMM_COLORPLOTS = 0;
  end
  if (mode == 1)
    EMHMM_COLORPLOTS = 1;
  end    
end

fprintf('** emhmm is set to %s figures\n', names{EMHMM_COLORPLOTS+1});