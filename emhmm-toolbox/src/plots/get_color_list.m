function [color] = get_color_list(numcolors)
% get_color_list - list of colors for each state
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-01-13
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong

% 2017-05-24: ABC - added input for number of colors, and added more colors > 7.
%                 - changed output to RGB colorspace, instead of letter codes.

% 2018-04-07: ABC - made yellow darker (from [1 1 0] to [0.85 0.85 0]) for visibility.
% 2018-05-22: v0.73 - made black into gray
% 2018-09-05: v0.74 - added option to use grayscale only

if nargin<1
  numcolors = 7;
end

global EMHMM_COLORPLOTS

if EMHMM_COLORPLOTS
  %% USE COLOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % if less than 7 colors, use existing color scheme.
  if numcolors <= 7
    %color = ['r','g','b','m', 'c', 'y', 'k'];
    %colorfull = {'red', 'green', 'blue', 'magenta', 'cyan', 'yellow', 'black'};
    
    color = {[1 0 0], [0 1 0], [0 0 1], [1 0 1], [0 1 1], [0.85 0.85 0], [0.3 0.3 0.3]};
    
  else
    
    % use jet colormap - flip it so that red is first.
    color = mat2cell(flipud(jet(numcolors)), ones(1,numcolors), 3);
  end

else 
  %% USE GRAYSCALE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  MAXGRAY = 1;
  MINGRAY = 0;
  
  color = mat2cell(linspace(MINGRAY, MAXGRAY, numcolors)'*[1 1 1], ones(1,numcolors), 3);

  % NOTE: Some gray levels may be hard to see from the underlying image.
  %       To override the gray cale levels, you can uncomment the below code
  %       and add/change gray levels as required.
  %color = {[0 0 0], [0.5 0.5 0.5], [1 1 1]};
  
end