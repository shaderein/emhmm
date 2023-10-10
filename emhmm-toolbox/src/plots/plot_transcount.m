function plot_transcount(M)
% plot_transcount - plot transition counts
%
%  plot_transcount(M)
%
% INPUT: 
%   M = transition count matrix (from vbhmm)
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-01-13
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong

% 2018-03-12: set min value to 0 for imagesc scaling
% 2018-5    : color text
% 2019-02-20: v0.75 - text scales with figure

maxM = max(M(:));
K = size(M,1);

% v0.75 - normalized font size
nfontsize = 0.25/max(3,K);
nfontsize2 = 0.25/3;

imagesc(M, [0 maxM]);
colorbar
set(gca, 'XTick', [1:K]);
set(gca, 'YTick', [1:K]);
for j=1:size(M,1)
  for k=1:size(M,2)
    mycolor = get_text_color(M(j,k) / maxM);
    text(k,j, sprintf('%.1f', M(j,k)), 'HorizontalAlignment', 'center', ...
      'FontUnits', 'normalized', 'FontSize', nfontsize, 'Color', mycolor);
  end
end
colormap gray
title('transition counts');
xlabel('to ROI');
ylabel('from ROI');
set(gca, 'FontUnits', 'normalized', 'FontSize', nfontsize2);
