function vhem_plot(group_hmms, imgfile, mode, subplotinfo, group_names)
% vhem_plot - plot group hmms from vhem_cluster
%
%   vhem_plot(group_hmms, imgfile, mode, subplotinfo, group_names)
%
% INPUTS
%   group_hmms = group_hmms output from vhem_cluster
%   imgfile    = filename of image for plotting under the fixations (optional)
%                or an actual image
%   mode       = 'c' -- use compact plots [default]
%              = 'b' -- use compact plot,  with transition matrix below
%              = 'h' -- HMM plots, each row shows one group HMM
%              = 'v' -- HMM plots, each column shows one group HMM
%  subplotinfo = use an existing figure and subplots for compact mode
%              = [py, px, ind1, ind2, ...]
%                 py = subplot y-size
%                 px = subplot x-size
%                 ind1, ind2, ... = subplot indices for each group
%              = [] - create new figure (default)
%  group_names = cell array of group names (default = {'Group 1', ...}
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-01-13
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong

% VERSIONS
%  2017-01-19 - updated for duration
%             - added 'h' and 'v' options for plotting
% 2017-08-11 - added subplotinfo option
% 2018-05-15 - v0.72 - add group_names
% 2018-06-14 - v0.74 - add support for passing an image directly
% 2020-03-21 - v0.78 - added 'b' option for compact-plot

if nargin<2
  imgfile = '';
end
if nargin<3
  mode = 'c';
end
if nargin<4
  subplotinfo = [];
end
if nargin<5
  group_names = {};
end

if ~isempty(imgfile)
  if ischar(imgfile)
    img = imread(imgfile);
  else
    img = imgfile;
  end
else
  img = [];
end

for i=1:length(group_hmms.hmms)
  if ~isempty(group_hmms.hmms{i})
    D = length(group_hmms.hmms{i}.pdf{1}.mean);
    break;
  end
end
if (D==2)
  nump = 3;
else
  nump = 4;
end

if ~isempty(subplotinfo) && ((mode ~= 'c') && (mode ~= 'b'))
  error('subplotinfo only supported with compact mode');
end

% plot size
G = length(group_hmms.hmms);
switch(mode) 
  case {'c', 'b'}
    % compact plot
    if isempty(subplotinfo)
      py = G;
      px = 1;
      plotselector = @(g,ind) g;
    else
      py = subplotinfo(1);
      px = subplotinfo(2);
      plotselector = @(g,ind) subplotinfo(g+2);
    end
    
  case 'h'
    % plot HMM as a row
    py = G;
    px = nump;
    plotselector = @(g,ind) ind+px*(g-1);
  case 'v'
    % plot HMM as a column
    py = nump;
    px = G;
    plotselector = @(g,ind) g+px*(ind-1);
end

if isempty(group_names)
  for g=1:G
    group_names{g} = sprintf('Group %d', g);
  end
end

if isempty(subplotinfo)
  figure
end

for g=1:G
  %title(['prior ' glab], 'interpreter', 'none');
  
  switch(mode)
    case {'c', 'b'}
      ind = 1;
      subplot(py,px,plotselector(g,ind))
      if ~isempty(group_hmms.hmms{g})
        if (mode == 'c')          
          vbhmm_plot_compact(group_hmms.hmms{g}, imgfile);
        else
          vbhmm_plot_compact(group_hmms.hmms{g}, imgfile, 'b');          
        end
      end
      if isfield(group_hmms, 'group_size')
        tmp = sprintf(' (size=%d)', group_hmms.group_size(g));
      else
        tmp = '';
      end
      title([group_names{g} tmp]);
      
    case {'v', 'h'}
      ind = 1;
      subplot(py,px,plotselector(g,ind))
      plot_emissions([], [], group_hmms.hmms{g}.pdf, img)
      if isfield(group_hmms, 'group_size')
        tmp = sprintf(' (size=%d)', group_hmms.group_size(g));
      else
        tmp = '';
      end
      title([group_names{g} tmp]);
      ind = ind+1;
      
      if (D>2)
        subplot(py,px,plotselector(g,ind))
        plot_emissions_dur([], [], group_hmms.hmms{g}.pdf);
        ind = ind+1;
      end
      
      subplot(py,px,plotselector(g,ind))
      plot_transprob(group_hmms.hmms{g}.trans)
      ind = ind+1;
    
      subplot(py,px,plotselector(g,ind))
      plot_prior(group_hmms.hmms{g}.prior)
      ind = ind+1;
  
    otherwise
      error('bad option');
  end
end
