function [] = savefigs(filename, formats, hfig)
% savefigs -- save the current figure, and export as an image
%
%   USAGE: savefigs(name, formats, hfig)
%
%   saves the current figure, as it appears on the screen,
%   and also exports as an EPS, PNG, or JPG file.
%     
%   INPUTS:
%      name - name to use for saving the file
%   formats - which formats to save as: fig, eps, png, jpg
%              default={'fig', 'png'}
%      hfig - the figure handle (default= the current figure)
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2019-07-13
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong


% 2018-07-17: add support for newer versions of MATLAB that preserve
%             figure size
% 2020-03-12: v0.77 - initial version

if nargin<2
  formats = {'fig', 'png'};
end

if exist('hfig', 'var')==0
    hfig = gcf;
end

tmp = version('-release');
matver = str2double(tmp(1:4));

% for MATLAB versions before 2016a: 
if matver < 2016
  % To set a figure to be printed the same size as shown on screen
  try
    set(hfig,'PaperPositionMode','auto');
  catch
  end
end

% save the figure
if any(strcmp('fig', formats))
  try
    saveas(hfig, [filename '.fig'], 'fig');
  catch
    warning('could not save figure as .fig. This is a bug in the Mac version of Matlab.');
  end
end

% save as EPS color file
if any(strcmp('eps', formats));
  saveas(hfig, [filename '.eps'], 'epsc');
end

% save as PNG (color) file
if any(strcmp('png', formats))
  saveas(hfig, [filename '.png'],'png');
end


% save as JPG (color) file
if any(strcmp('jpg', formats))
  saveas(hfig, [filename '.jpg'],'jpg');
end


