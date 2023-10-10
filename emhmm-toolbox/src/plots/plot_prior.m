function plot_prior(prior, dp)
% plot_prior - plot fixations on an image
%
%  plot_prior(prior)
%
% INPUT: 
%   prior - prior distribution
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-01-13
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong

% 2017-05-24: ABC - updated to use new color_list
% 2019-02-20: v0.75 - text changes with figure size
% 2019-06-29: v0.77 - show prior probabilites
if (nargin<2)
  dp = 2;
end

K = length(prior);

nfontsize2 = 0.25/3;

color = get_color_list(K);

set(gca, 'XTick', [1:K]);
hold on
for i=1:K
  bar(i, prior(i), 'FaceColor', color{i});  
  
  if (prior(i)<0.5) vertalign = 'bottom';
  else vertalign = 'top';
  end
  text(i,prior(i), prob2str(prior(i), dp), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', vertalign, ...
    'FontUnits', 'normalized', 'FontSize', nfontsize2);
end
grid on
title('prior');
ylabel('probability');
xlabel('ROI');
set(gca, 'FontUnits', 'normalized', 'FontSize', nfontsize2);
