function [hfig, roiseq] = vhem_plot_fixations(data, group_hmms, imgfile, opt, subsample, figopts)
% vhem_plot_fixations - plot fixations for group_hmms
%
%   vhem_plot_fixations(data, group_hmms, imgfile, opt, subsample)
%
% INPUTS
%   data       = data cell (used to learn individual HMMs)
%   group_hmms = group_hmms output from vhem_cluster
%   imgfile    = filename of image for plotting under the fixations (optional)
%                or a loaded image.
%   opt        = 'g' - show fixations for each group [default]
%              = 's' - separately show fixations from each subject
%              = 'c' - use compact plot to show ROI/fixations and transition matrix
%              = 'b' - use compact plot w/ transition matrix on bottom
%    subsample = sub-sample the data when plotting (default=1, no subsampling)
%                 2 = plot 1/2 of the data 
%                 3 = plot 1/3 of the data, etc
%                -N = plot N total fixation points
% OUTPUTS
%       hfig = handle to figure
%  roiseq{i} = cell array of state sequences for data{i} under 
%              its most likely group HMM (i.e., group_hmms.label(i))
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-08-22
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong

% VERSIONS
% 2018-06-14: v0.74 - added output fig handle, support for passing image directly
% 2018-11-24: v0.74 - handle empty data, fix bug when concat data for viewing
% 2019-01-22: v0.75 - added state sequence output
% 2019-01-22: v0.75 - added compact plot option
% 2019-02-20: v0.75 - better margins and layout, plot duration data
% 2019-02-21: v0.75 - subsample fixations
% 2020-03-22: v0.78 - added 'b' plot
% 2021-07-30: v0.80 - allow different S for each group

if nargin<3
  imgfile = '';
end
if nargin<4
  opt = 'g';
end
if nargin<5
  subsample = 1;
end
if nargin<6
  figopts = {};
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

K = length(group_hmms.hmms);
D = length(group_hmms.hmms{1}.pdf{1}.mean);

gamma = cell(1,length(data));

if (nargout >= 2)
  roiseq = cell(1,length(data));
end

% for each subject
for i=1:length(data)
  mydata = data{i};  
  
  % 2018-11-24: v0.74 - handle empty data
  if isempty(mydata)
    mygamma = {};
    ss = {};
  else
    
    k = group_hmms.label(i);
    
    % v0.80 - allow S to be different for each cluster
    S = length(group_hmms.hmms{k}.prior);
    
    % get most-likely ROI sequences
    myhmm = group_hmms.hmms{k};
    ss = vbhmm_map_state(myhmm, mydata);
    
    % setup as gamma for plot_emissions
    mygamma = cell(1,length(mydata));
    for j=1:length(mydata)
      mygamma{j} = zeros(S, size(mydata{j},1));
      for q=1:size(mydata{j},1)
        mygamma{j}(ss{j}(q), q) = 1;
      end
    end
  end
  
  gamma{i} = mygamma; %TODO: classified cluster for each fixation point?
  
  if nargout >= 2
    roiseq{i} = ss;
  end
end

hfig = figure(figopts{:});

% plot fixations for each group
if any(opt=='g') || any(opt=='c') || any(opt=='b')
  if (D==2)
    if any(opt=='g')
      py = 1;
      px = K;
      gap = [];
      margh = [];
      margw = [];
    else
      py = K;
      px = 1;
      gap = [0.04];
      margh = [];
      margw = [];
    end
    plotselector = @(y,x) y;    
  elseif (D==3)
    if any(opt=='g')
      py = 2;
      px = K;
      gap = [0.1 0.1];      
      margh = [0.1 0.05];
      margw = [0.1 0.05];
      plotselector = @(y,x) (x-1)*px + y;
    else
      py = K;
      px = 1;
      gap = [0.04];
      margh = [];
      margw = [];
      plotselector = @(y,x) y;
    end
  else
    error('dim>3 not supported yet');
  end
  
  for k=1:K
    myg = group_hmms.groups{k};
    % v0.74 - BUG FIX: cat in the correct direction
    if any(cellfun(@isrow, data))
      gdata = cat(2,data{myg});  % data is row cell vector
    else
      gdata = cat(1,data{myg}); % data is column cell vector
    end
    ggamma = cat(2,gamma{myg});
    
    if any(opt=='g')
      subtightplot(py,px,plotselector(k,1), gap, margh, margw)
      plot_emissions(gdata, ggamma, group_hmms.hmms{k}.pdf, img, 't', subsample);
      title(sprintf('Group %d', k))
      
      if (D==3)
        subtightplot(py,px,plotselector(k,2), gap, margh, margw)
        plot_emissions_dur(gdata, ggamma, group_hmms.hmms{k}.pdf, subsample);
      end
      
    elseif any(opt=='c') || any(opt=='b')
      subtightplot(py,px,plotselector(k,1), gap, margh, margw)
      tmphmm = group_hmms.hmms{k};
      tmphmm.gamma = ggamma;
      if any(opt=='c')
        mycode = 'r';
      else
        mycode = 'b';
      end
      vbhmm_plot_compact(tmphmm, img, mycode, gdata, subsample);
      title(sprintf('Group %d', k));
      
    else
      error('not sure how you got here.');
    end
      
  end
  
% plot fixations for each subject in a group
elseif(any(opt=='s'))
  if (D==2)
    py = ceil(length(data)/5);    
    px = 5;
    plotselector = @(y,x) y;
  elseif (D==3)
    py = ceil(2*length(data)/8);
    px = 8;    
    plotselector = @(y,x) (y-1)*2 + x;
  else 
    error('dim>3 not supported yet');
  end
  
  for i=1:length(data)
    gdata = data{i};
    ggamma = gamma{i};
    myg = group_hmms.label(i);
    
    subplot(py,px,plotselector(i,1))
    plot_emissions(gdata, ggamma, group_hmms.hmms{myg}.pdf, img, '', subsample);
    title(sprintf('Subj %d (Grp %d)', i, myg));
    
    if (D==3)
      subplot(py,px,plotselector(i,2))
      plot_emissions_dur(gdata, ggamma, group_hmms.hmms{k}.pdf, subsample);
    end
  end
end

