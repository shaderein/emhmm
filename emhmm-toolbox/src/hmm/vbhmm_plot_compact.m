function vbhmm_plot_compact(hmm, imgfile, plotmode, data, subsample)
% vbhmm_plot_compact - compact plot of HMM
%
%   vbhmm_plot_compact(hmm, imgfile, plotmode, data, subsample)
%
% INPUTS
%   hmm      = hmm from vbhmm_learn
%   imgfile  = filename of image for plotting under the fixations (optional)
%              or an actual image matrix.
%   plotmode = 'r' - plot the transition matrix to the right of ROI plot (default)
%            = 'b' - plot the transition matrix below the ROI plot
%               (if 'r' or 'b' are not specified, then no transition matrix is plotted)
%   data     = also plot the fixations (optional)
%  subsample = subsample the data before plotting (default=1)
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-01-13
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong

% 2016-11-22: ABC - created
% 2017-05-24: ABC - updated to use new color_list
% 2018-02-01: ABC - added option to plot fixations
% 2018-03-11: v0.72 - added support for groups
% 2018-06-13: v0.74 - can pass an image instead of imgfile.
% 2019-01-23: v0.75 - text font and line width now adapts to plot size
% 2019-02-22: v0.75 - more grayscale support, subsampling
% 2020-03-21: v0.78 - added 'b' option to show transition matrix below
%                   - make padding adapt to transition matrix size
%                   - make fontsize better for odd-shaped transition matrices.
%                   - only make grayscale colormap for local axis
% 2021-07-27: v0.80 - support empty image

global EMHMM_COLORPLOTS

if nargin<2
  imgfile = '';
end
if nargin<3
  plotmode = 'r';
end
if nargin<4
  data = [];
end
if nargin<5
  subsample = 1;
end

if ~isempty(imgfile)
  if ischar(imgfile)
    % image filename 
    img = imread(imgfile);
  else
    % an actual image
    img = imgfile;
  end
  
  if ndims(img)==2
    img = cat(3,img,img,img);
  end
else
  img = [];
end

% get fixation plot dimensions
if ~isempty(img)
  xmin_img = 0;
  xmax_img = size(img,2);
  ymin_img = 0;
  ymax_img = size(img,1);
else
  % v0.80 - support empty image
  tmpxs = zeros(0,2);
  for k=1:length(hmm.pdf)
    tmpm = hmm.pdf{k}.mean;
    tmps = sqrt(diag(hmm.pdf{k}.cov)');
    tmpxs(end+1,:) = tmpm-2.5*tmps;
    tmpxs(end+1,:) = tmpm+2.5*tmps;
  end
  xmin_img = min(tmpxs(:,1));
  xmax_img = max(tmpxs(:,1));
  ymin_img = min(tmpxs(:,2));
  ymax_img = max(tmpxs(:,2));
end

nfontsize = 0.25/4;

D = length(hmm.pdf{1}.mean);
K = length(hmm.pdf);

[colors] = get_color_list(K);

% get location of transition matrix
if any(plotmode == 'r')
  xmin_trans = xmax_img;
  xmax_trans = xmax_img + min(xmax_img, ymax_img);
  ymin_trans = ymin_img;
  ymax_trans = ymax_img;
  if (D>2)
    % set the location of the duration emissions
    xmin_dur = xmin_trans;
    xmax_dur = xmax_trans;
    ymin_dur = ymin_trans;
    ymax_dur = ymax_trans;
    
    % and move the transition matrix to the right
    w = (xmax_dur-xmin_dur);
    xmin_trans = xmin_trans + w;
    xmax_trans = xmax_trans + w;
  end
end

if any(plotmode == 'b')
  xmin_trans = 0;
  xmax_trans = xmax_img;
  ymin_trans = ymax_img;
  ymax_trans = ymax_img + (xmax_img - xmin_img)/2;
  
  if (D>2)
    error('not supported yet');
  end
end

xmax_all = max(xmax_img, xmax_trans);
ymax_all = max(ymax_img, ymax_trans);

% v0.78 - adapt padding
tmp = [ (xmax_trans - xmin_trans) / (K+0.5), ...
        (ymax_trans - ymin_trans) / (K+1.5) ];
allpadding = min(tmp);


%% plot the ROIs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(data)
  plot_emissions([], [], hmm.pdf, img, 'b', 1, [xmax_all, ymax_all]);
else
  plot_emissions(data, hmm.gamma, hmm.pdf, img, 'b', subsample, [xmax_all, ymax_all]);
end
title('');

if isempty(plotmode)
  return
end
hold on

%% plot the duration emissions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (D>2)
  if EMHMM_COLORPLOTS
    lw = 2;
    bgcol = 'w';
  else
    lw = 4;
    bgcol = 0.5*(colors{end}+colors{end-1});
  end

  
  % get size
  xw = xmax_dur - xmin_dur;
  yh = ymax_dur - ymin_dur;
  
  % get line for each pdf
  dur_z    = {};
  dur_t    = {};
  dur_mu_z = {};
  legs = {};
  for k=1:K
    mu     = hmm.pdf{k}.mean(3);
    sigma2 = hmm.pdf{k}.cov(3,3);
    ss = sqrt(sigma2);
    tmin = max(mu - 3*ss, 0);
    tmax = mu + 3*ss;
    t = linspace(tmin, tmax, 100);
    z = normpdf(t, mu, ss);
    
    dur_t{k}    = t;
    dur_z{k}    = z;
    dur_mu_z{k} = [mu, normpdf(mu,mu,ss)];
    
    if EMHMM_COLORPLOTS
      legs{k} = sprintf('\\color[rgb]{%0.3f,%0.3f,%0.3f}%d: %d\\pm%d', ...
        colors{k}(1), colors{k}(2), colors{k}(3), k, round(mu), round(ss));
    else
      legs{k} = sprintf('%d: %d\\pm%d', k, round(mu), round(sqrt(sigma2)));
    end
    
  end

  % find the t range and pt range
  tmin = 0;
  tmax = max(cat(2, dur_t{:}));
  zmin = 0;
  zmax = max(cat(2, dur_z{:}));

  % remap plot to canvas (image) coordinates
  padding = allpadding / 2;
  yaxpadding = 65;
  textpadding = 5;
  textoffset = 20;
  t_map = @(t) ((t-tmin) / (tmax-tmin))*(xw-padding*2) + xmin_dur+padding;
  z_map = @(z) ymax_dur-yaxpadding - ((z-zmin) / (zmax-zmin))*(yh-yaxpadding);  % upside down
  
  % plot axes
  fill([xmin_dur+padding, xmin_dur+padding, xmax_dur-padding, xmax_dur-padding], ...
       [ymin_dur, ymax_dur-yaxpadding, ymax_dur-yaxpadding, ymin_dur], ...
       bgcol, 'linewidth', 0.5);
  
  % v0.75 fill bottom if data
  if ~isempty(data)    
    fill([xmin_dur+padding, xmin_dur+padding, xmax_dur-padding, xmax_dur-padding], ...
      [ymax_dur-yaxpadding, ymax_dur, ymax_dur, ymax_dur-yaxpadding], ...
      bgcol, 'linewidth', 0.5);
  end
    
  % plot axes text values
  ytext = ymax_dur-yaxpadding+textpadding;
  text(xmin_dur+padding, ytext, sprintf('%d', floor(tmin)), ...
    'HorizontalAlignment', 'left', 'FontUnits', 'normalized', 'FontSize', 0.75*nfontsize, ...
    'VerticalAlignment', 'top');
  text(xmax_dur-padding, ytext, sprintf('%d', ceil(tmax)), ...
    'HorizontalAlignment', 'right', 'FontUnits', 'normalized', 'FontSize', 0.75*nfontsize, ...
    'VerticalAlignment', 'top');  
  
  
  % v0.75 - plot duration data
  if ~isempty(data)    
    out = cat(1, data{:});
    gamma = cat(2, hmm.gamma{:});
    [~,cl] = max(gamma,[],1);
          
    % v0.75 - subsampling
    if subsample < 0
      % plot N values
      ss = floor(length(cl) / (-subsample));
    elseif subsample > 0
      % fixed subsampling
      ss = subsample;
    else
      % no subsampling
      ss = 0;
    end
    
    for k=1:K
      if EMHMM_COLORPLOTS
        mycol = colors{k};
        scatteropts = {};
        mysize = 6;
      else
        mycol = 'k';
        scatteropts = {'MarkerFaceColor', colors{k}};
        mysize = 20;
      end
      
      mytext = ytext+(k-1)*textoffset;
      yrand  = textoffset*0.8;
      
      ii = find(cl==k);           
      % v0.75 - subsampling
      if (ss ~= 0)
        ii = ii(1:ss:end);
      end
      
      % randomly y-offset the data
      randoffset = yrand*(rand(1,length(ii)));
      scatter(t_map(out(ii,3)), mytext*ones(1,length(ii)) + randoffset, mysize, mycol, 'o', scatteropts{:});
    end
  end
  
  
  % plot pdfs
  for k=1:K
    mytext = ytext+(k-1)*textoffset;
    plot(t_map(dur_t{k}), z_map(dur_z{k}),'Color', colors{k}, 'linewidth', lw);
    
    % plot center line
    plot(t_map(dur_mu_z{k}(1))*[1 1], [mytext, z_map([dur_mu_z{k}(2)])], '--', 'Color', colors{k}, 'linewidth', lw/2);
    
    if isempty(data)
      % plot duration
      text(t_map(dur_mu_z{k}(1)), mytext, sprintf('%d', round(dur_mu_z{k}(1))), ...
        'HorizontalAlignment', 'center', 'FontUnits', 'normalized', 'FontSize', 0.75*nfontsize, ...
        'VerticalAlignment', 'top');
    end
  end
  
  % plot labels
  for k=1:K     
    if EMHMM_COLORPLOTS
      textcol = colors{k};
      textopts = {};
    else
      % v0.75 grayscale
      textopts = {'FontWeight', 'Bold', 'BackgroundColor', colors{k}};
      textcol = get_text_color(colors{k});
    end
    
    text(t_map(dur_mu_z{k}(1)), z_map(dur_mu_z{k}(2)/2), sprintf('%d', k), 'color', textcol, ...
      'horizontalalignment', 'center', 'FontUnits', 'normalized', 'FontSize', nfontsize, textopts{:});
  end
  
  % more compact legend
  if EMHMM_COLORPLOTS
    text(xmax_dur-padding, ymin_dur, legs, 'color', 'black', ...
      'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
      'BackgroundColor', 'white', ...
      'FontUnits', 'normalized', 'FontSize', nfontsize*0.75, ...
      'EdgeColor', 'black', 'Margin', 1);
    
  else
    % v0.75 grayscale mode
    ytmp = ymin_dur;
    for k=1:K
      H = text(xmax_dur-padding, ytmp, legs{k}, 'color', get_text_color(colors{k}), ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
        'BackgroundColor', colors{k}, ...
        'FontUnits', 'normalized', ...
        'FontSize', nfontsize*0.75, 'EdgeColor', 'black', 'Margin', 1);
      ytmp = ytmp+H.Extent(4);
    end
  end
  
  
end

if ~iscell(hmm.prior)
  mypriors = {hmm.prior};
  mytranss = {hmm.trans};
else
  % groups
  mypriors = hmm.prior;
  mytranss = hmm.trans;
end

%padding=40;
padding = allpadding;

  
%nfontsize = 0.8 * ((xmax_trans-xmin_trans) / max(3,K) ) / xmax_trans;


for g = 1:length(mypriors)
  myprior = mypriors{g};
  mytrans = mytranss{g};
  if (g>1)
    diff = xmax_trans - xmin_trans;
    xmin_trans = xmax_trans;
    xmax_trans = xmax_trans + diff;
  end
  
  %% plot the prior & transition matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if any(plotmode == 'r')
    xpadding = padding;
    ypadding = 0;
  end
  if any(plotmode == 'b')
    xpadding = padding;
    ypadding = padding/4;
  end
  
  % x spacing (common for both)
  tx = linspace(xmin_trans+xpadding, xmax_trans, K+1);
  tx = 0.5*(tx(2:end)+tx(1:end-1));
  
  % overall y spacing
  ty = linspace(ymin_trans+padding+ypadding, ymax_trans, K+2);
  ty = 0.5*(ty(2:end)+ty(1:end-1));
  
  % y spacing for prior
  typ = ty(1)-padding;
  
  % y spacing for transition matrix
  tyt = ty(2:end);
  
  % cell size
  if length(tyt)>1
    dy = (tyt(2) - tyt(1));
  else
    dy = (ymax_trans-ymin_trans-padding);
  end
  if length(tx)>1
    dx = (tx(2) - tx(1));
  else
    dx = (xmax_trans-xmin_trans-padding);
  end
  
  % v0.78 - better fontsize
  tmp = [(dy / ymax_all), ...   % using cell height
         (dx/1.5) / ymax_all ];  % using cell width (.XX is 1.5 times the font height).
  
  nfontsize = min(tmp);
  
  % plot prior (hack to trick matlab into plotting a vector)
  %  ... make it plot a matrix instead.
  imagesc(tx, typ+[-dy/4, dy/4], [myprior(:)'; myprior(:)'], [0 1]);
  
  imagesc(tx, tyt, mytrans, [0 1]);

  %colorbar
  
  % plot probabilities
  for j=1:K
    mycolor = get_text_color(myprior(j));
    text(tx(j),typ(1), prob2str(myprior(j), 2), ...
      'HorizontalAlignment', 'center', 'FontUnits', 'normalized', ...
      'FontSize', nfontsize, 'Color', mycolor);
  end
  
  for j=1:K
    for k=1:K
      mycolor = get_text_color(mytrans(j,k));
      text(tx(k),tyt(j), prob2str(mytrans(j,k), 2), ...
        'HorizontalAlignment', 'center', 'FontUnits', 'normalized', ...
        'FontSize', nfontsize, 'Color', mycolor);
    end
  end
  
  %nfontsize = 3.5*(padding*3/4) / xmax_trans;
  %nfontsize = 2*(padding*3/4) / xmax_trans;

  % v0.78 - adapt fontsize better
  nfontsize = 0.8*(3/4*padding) / ymax_all;  % for rows
  nfontsizec = min(nfontsize, 0.8*(dx/2) / ymax_all);  % for columns
  nfontsizep = 0.8*(3/4*padding / 1.5) / ymax_all;  % for prior
  
  % plot ROI color strips
  for j=1:K
    mycolor = get_text_color(colors{j});
        
    % labels for rows
    rectangle('Position', [xmin_trans+padding/4, tyt(j)-dy/2, padding*3/4, dy], 'FaceColor', colors{j});
    text(xmin_trans+padding*2.5/4, tyt(j), sprintf('%d', j), ...
      'HorizontalAlignment', 'center', 'FontUnits', 'normalized', 'FontSize', nfontsize, ...
      'Color', mycolor);
    
    % labels for columns
    rectangle('Position', [tx(j)-dx/2, tyt(1)-dy/2-padding*3/4, dx, padding*3/4], 'FaceColor', colors{j});
    text(tx(j), tyt(1)-dy/2-padding*1.5/4, sprintf('to %d', j), ...
      'HorizontalAlignment', 'center', 'FontUnits', 'normalized', 'FontSize', nfontsizec, ...
      'Color', mycolor);
  end
  
  % prior
  text(xmin_trans+padding*2.5/4, typ, 'prior', ...
    'Rotation', 90, 'FontUnits', 'normalized', 'FontSize', nfontsizep, 'HorizontalAlignment', 'center')
  
end


% reset axis
axis([min(xmin_img, xmin_trans), ...
      max(xmax_img, xmax_trans), ...
      min(ymin_img, ymin_trans), ...
      max(ymax_img, ymax_trans)]);
 
% v0.78 - only change colormap of local axis
colormap(gca, gray)
%colorbar;

%plot_transprob(hmm.trans)
%plot_prior(hmm.prior)

