function plot_emissions_dur(data, gamma, pdf, subsample)
% plot_emissions_dur - plot emission densities for durations
%
%     plot_emissions_dur(data, gamma, pdf)
%
%  INPUTS
%    data  = data (same format as used by vbhmm_learn)
%    gamma = assignment variables for each data point to emission density
%    pdf   = the emission densities (from vbhmm)
%   subsample = subsampling of data for visualization
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-01-13
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong


% VERSIONS
% 2017-01-17: ABC - iniital version
% 2017-05-24: ABC - updated to use new color_list
% 2019-02-20: v0.75 - text changes with figure size
% 2019-02-22: v0.75 - grayscale support, subsampling

global EMHMM_COLORPLOTS

if nargin<4
  subsample = 1;
end

hold on

nfontsize = 0.25/3;
nfontsize2 = 0.25/3;

maxzs = [];
maxt = 0;
K = length(pdf);
color = get_color_list(K);

if EMHMM_COLORPLOTS
  lw = 2;
else
  lw = 4;
  set(gca, 'Color', 0.5*(color{end}+color{end-1}));
end

legs = {};
for k = 1:K
  mu    = pdf{k}.mean(3);
  sigma2 = pdf{k}.cov(3,3);
  [t, z] = plot1D(mu, sigma2, color{k}, lw);
  maxzs(k) = max(z);
  maxt = max(max(t), maxt);
  if EMHMM_COLORPLOTS
    legs{k} = sprintf('\\color[rgb]{%0.3f,%0.3f,%0.3f}%d: %d\\pm%d', ...
      color{k}(1), color{k}(2), color{k}(3), k, round(mu), round(sqrt(sigma2)));
  else
    legs{k} = sprintf('%d: %d\\pm%d', k, round(mu), round(sqrt(sigma2)));
  end
end
maxz = max(maxzs);

% plot data
if ~isempty(data)      
  out = cat(1, data{:});
  gamma = cat(2, gamma{:});
  [~,cl] = max(gamma,[],1);
  K = size(gamma,1);    
  
  D = size(data{1}, 2);
  if D <= 2
    error('data has no duration (dimension 3)');    
  end
  
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
  
  yoffset = maxz/10;
  yrand   = yoffset/4;
  for k=1:K
    if EMHMM_COLORPLOTS
      mycol = color{k};
      scatteropts = {};
      mysize = 6;
    else
      mycol = 'k';
      scatteropts = {'MarkerFaceColor', color{k}};
      mysize = 20;
    end
    
    ii = find(cl==k);
    
    % v0.75 - subsampling
    if (ss ~= 0)
      ii = ii(1:ss:end);
    end  
    
    % randomly y-offset the data
    randoffset = yrand*(2*rand(1,length(ii))-1);
    scatter(out(ii,3), -yoffset*(k-1)*ones(1,length(ii)) + randoffset, mysize, mycol, 'o', scatteropts{:});
  end    
  
  ymin = -yoffset*((k-1)+0.5);
else
  ymin = 0;
end

% plot labels
for k = 1:K
  if EMHMM_COLORPLOTS
    textcol = color{k};
    textopts = {};    
  else
    % v0.75 grayscale
    textopts = {'FontWeight', 'Bold', 'BackgroundColor', color{k}};
    textcol = get_text_color(color{k});   
  end
  
  mu = pdf{k}.mean(3);
  text(mu, maxzs(k)/2, sprintf('%d', k), 'color', textcol, ...
    'horizontalalignment', 'center', ...
    'FontUnits', 'normalized', 'FontSize', nfontsize, textopts{:});
end

% get axes
ymax = maxz;
xmin = 0;
xmax = maxt;

% plot axes
plot([xmin, xmax], [0, 0], 'k-');

% compact legend
if EMHMM_COLORPLOTS
  text(xmax, ymax, legs, 'color', 'black', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'BackgroundColor', 'white', ...
    'FontUnits', 'normalized', ...
    'FontSize', nfontsize*0.75, 'EdgeColor', 'black', 'Margin', 1);
  
else
  % v0.75 grayscale mode
  ymaxtmp = ymax;
  for k=1:K
    H = text(xmax, ymaxtmp, legs{k}, 'color', get_text_color(color{k}), ...
      'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
      'BackgroundColor', color{k}, ...
      'FontUnits', 'normalized', ...
      'FontSize', nfontsize*0.75, 'EdgeColor', 'black', 'Margin', 1);
    ymaxtmp = ymaxtmp-H.Extent(4);
  end
end


grid on
hold off
title('ROIs (duration)');
xlabel('t');
ylabel('p(t)');

% reset axes
axis([xmin, xmax, ymin, ymax]);

set(gca, 'FontUnits', 'normalized', 'FontSize', nfontsize2);



%% plot a Gaussian curve
function [t, z, h] = plot1D(mu, sigma2, color, lw)

ss = sqrt(sigma2);

tmin = mu - 3*ss;
tmax = mu + 3*ss;
t = linspace(tmin, tmax, 100);
z = normpdf(t, mu, ss);

h = plot(t, z,'Color', color, 'linewidth', lw);
plot([mu, mu], [0, normpdf(mu,mu,ss)], '--', 'Color', color, 'linewidth', lw/2);
