function plot_emissions(data, gamma, pdf, img, opt, subsample, plotsize)
% plot_emissions - plot emission densities
%
%  Plots ellipses for each ROI. The ellipse is 2 standard-deviations (95% percentile).
%
%     plot_emissions(data, gamma, pdf, img, opt, subsample, plotsize)
%
%  INPUTS
%    data  = data (same format as used by vbhmm_learn)
%    gamma = assignment variables for each data point to emission density
%    pdf   = the emission densities (from vbhmm)
%    img   = image for background
%    opt   = plotting options:
%            's' - show sequence lines
%            't' - use transparency with markers
%            'b' - use thicker lines, and text scales with image
%    subsample = sub-sample the data when plotting (default=1, no subsampling)
%                 2 = plot 1/2 of the data 
%                 3 = plot 1/3 of the data, etc
%                -N = plot N total fixation points
%               not valid when using 's' option
%    plotsize = the [x,y] size of the plot (used for making adaptive text) (optional)
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-01-13
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong


% VERSIONS
% 2017-01-17: ABC - add duration as marker size.
% 2017-05-24: ABC - updated to use new color_list
% 2017-06-21: ABC - use our own "showim" instead of "imshow"
% 2017-08-22: ABC - added option to show sequence
% 2018-02-01: ABC - added warning if no transparency
% 2018-05-16: v0.72 - remove dependence on stats toolbox
% 2018-09-05: v0.74 - grayscale plot: increase ROI thickness, bold text
% 2018-10-12: v0.74 - use MarkerEdgeAlpha; added sub-sampling option.
% 2019-02-20: v0.75 - text size changes with figure size
% 2020-03-21: v0.78 - fixed text size for wide images, added plotsize

global EMHMM_COLORPLOTS

if (nargin<5)
  opt = '';
end
if nargin<6
  subsample = 1;
end
if nargin<7
  plotsize = [];
end

%nfontsize2 = 0.25/3;
color = get_color_list(length(pdf));

if ~isempty(img)
  showim(img);
end
axis ij
hold on

if ~isempty(data)
  out = cat(1, data{:});
  gamma = cat(2, gamma{:});
  [~,cl] = max(gamma,[],1);
  K = size(gamma,1);
  
  D = size(data{1}, 2);
  if D > 2
    myscale = get_duration_rescale(data);
  end

  if any(opt=='s')
    if subsample ~= 1
      error('cannot subsample when showing sequence lines')
    end
    for n=1:length(data)
      plot(data{n}(:,1), data{n}(:,2), '-k');
    end
  end
  
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
    ii = find(cl==k);
    if (ss ~= 0)
      ii = ii(1:ss:end);
    end
            
    if (D == 2)
      h = scatter(out(ii,1),out(ii,2), 6, color{k}, 'o');
    else      
      h = scatter(out(ii,1),out(ii,2), myscale(out(ii,3)),  ...
        color{k}, 'o');
    end
    
    if any(opt=='t')
      if isprop(h, 'MarkerEdgeAlpha')
        h.MarkerEdgeAlpha = 0.7;
      end
      if isprop(h, 'MarkerFaceAlpha')
        h.MarkerFaceAlpha = 0.7;
      else
        % in MATLAB < 2015, alpha markers have not been implemented yet
        % so hack it.
        delete(h)
        % replot
        if (D == 2)
          h = scatter2(out(ii,1),out(ii,2), 6, color{k}, 'o');
        else
          h = scatter2(out(ii,1),out(ii,2), myscale(out(ii,3)),  ...
            color{k}, 'o');
        end
        warning('plot_emissions:transparency', 'Your MATLAB version cannot use transparency for the fixation markers. The plots will not look nice.');
        % as a courtesy, turn off this warning now
        warning('off', 'plot_emissions:transparency');
      end
    end
    
    
  end
else
  K = length(pdf);
end

for k = 1:K
  mu = pdf{k}.mean;
  sigma(:,:) = pdf{k}.cov;
  plot2D(mu, sigma, color{k}, opt)
  %text(mu(1), mu(2), sprintf('%d', k));
end

% v0.78 - adapt text size better for long images
nfontsize = 0.1; % default
if ~isempty(img)
  tmp = [0.05*size(img,2), 0.1*size(img,1)];
  if ~isempty(plotsize)
    sc = plotsize(2);
  else
    sc = size(img,1);
  end
  nfontsize = max(tmp) / sc;
end

for k = 1:K
  if EMHMM_COLORPLOTS
    textopts = {};
    textcol = color{k};
    if any(opt=='b')
      textopts = {'FontWeight', 'Bold', ... %'BackgroundColor', 'w', ...
        'FontUnits', 'Normalized', 'FontSize', nfontsize};
    end
    
  else
    textopts = {'FontWeight', 'Bold', 'BackgroundColor', color{k}, ...
      'FontUnits', 'Normalized', 'FontSize', nfontsize};
    textcol = get_text_color(color{k});    
  end
 
  mu = pdf{k}.mean;
  
  if any(opt=='b')
    % v0.75 - bold text -- make this its own function
    offsets = 2*[[-1, 1]; [0, 1]; [1, 1]; ...
                [-1, 0]; [1, 0]; ...
                [-1, -1]; [0, -1]; [1, -1]];
    
    for ii = 1:size(offsets,1)
      text(mu(1)+offsets(ii,1), mu(2)+offsets(ii,2), sprintf('%d', k), 'color', 'w', ...
        'horizontalalignment', 'center', textopts{:});
    end
  end
  
  text(mu(1), mu(2), sprintf('%d', k), 'color', textcol, ...
    'horizontalalignment', 'center', textopts{:});
  
  
end
hold off
title('ROIs');
%set(gca, 'FontUnits', 'normalized', 'FontSize', nfontsize2);


%% plot a Gaussian as ellipse
function plot2D(mu, Sigma, color, opt)

global EMHMM_COLORPLOTS
if EMHMM_COLORPLOTS
  lw = 1; % color
  if any(opt=='b')
    lw = 2;
  end
else
  lw = 2; % grayscale
end


% truncate to 2D
mu = mu(1:2);
Sigma = Sigma(1:2,1:2);

mu = mu(:);
if ~any(isnan(Sigma(:))) && ~any(isinf(Sigma(:)))
  [U,D] = eig(Sigma);
  n = 100;
  t = linspace(0,2*pi,n);
  xy = [cos(t);sin(t)];
  %k = sqrt(conf2mahal(0.95,2));
  k = sqrt(5.9915);
  w = (k*U*sqrt(D))*xy;
  z = repmat(mu,[1 n])+w;
  h = plot(z(1,:),z(2,:),'Color',color,'LineWidth',lw);
end

%function m = conf2mahal(c,d)
%m = chi2inv(c,d);

function h = scatter2(X, Y, S, C, m)

for i=1:length(X)
  cx = X(i);
  cy = Y(i);
  if length(S) ~= 1
    cs = S(i)/4;
  else
    cs = S/4;
  end
  
  patch([cx-cs, cx+cs, cx+cs, cx-cs, cx-cs], [cy-cs,cy-cs,cy+cs,cy+cs,cy-cs], ...
    C, 'FaceAlpha', 0.7, 'EdgeColor', 'none');
end
 h = [];