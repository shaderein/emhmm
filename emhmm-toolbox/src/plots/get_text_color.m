function mycolor = get_text_color(rgb)
% get_text_color - get the color for text with rgb background.
%
%  rgb = a RGB color vector [r, g, b].
%      = or a gray level between 0 and 1.
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-01-13
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong

if length(rgb) == 1
  p = rgb;
else
  %p = rgb2gray(rgb);
  %p = p(1);
  p = 0.2989 * rgb(1) + 0.5870 * rgb(2) + 0.1140 * rgb(3);
end
  
if (p > 0.3)
  mycolor = 'k';
else
  mycolor = 'w';
end
