function plot_transprob(trans, dp)
% plot_transprob - plot transition matrix
%
%  plot_transprob(M, dp)
%
% INPUT: 
%   M = transition matrix (from vbhmm)
%  dp = number of decimal places [default=2]
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-01-13
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong

% 2017-05-24: ABC - remove leading 0 from probabilities
% 2018-02-13: ABC - color-coded state labels
% 2019-02-20: v0.75 - text scales with figure

if (nargin<2)
  dp = 2;
end

K = size(trans,1);

% v0.75 - normalized font size
nfontsize = 0.33/max(3,K);
nfontsize2 = 0.25/3;

imagesc(trans, [0 1])

set(gca, 'XTick', [1:K]);
set(gca, 'YTick', [1:K]);
for j=1:size(trans,1)
  for k=1:size(trans,2)
    mycolor = get_text_color(trans(j,k));
    text(k,j, prob2str(trans(j,k), dp), ...
      'HorizontalAlignment', 'center', ...
      'FontUnits', 'normalized', 'FontSize', nfontsize, ...
      'Color', mycolor);
  end
end

% color tick labels
[color] = get_color_list(K);
TickLabels = {};
for j=1:K
  TickLabels{j} = sprintf('\\bf\\color[rgb]{%g,%g,%g}%d', ...
    color{j}(1), color{j}(2), color{j}(3), j);
end
set(gca, 'XTickLabel', TickLabels);
set(gca, 'YTickLabel', TickLabels);
colormap(gca, gray)
title('transition matrix');
xlabel('to ROI');
ylabel('from ROI');

set(gca, 'FontUnits', 'normalized', 'FontSize', nfontsize2);

colorbar
