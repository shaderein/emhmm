function plot_emcounts(N)
% plot_emcounts - plot emission counts
%
%  plot_emcounts(N) 
%
% INPUTS
%   N = emmission histogram
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-01-13
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong

% 2017-05-24: ABC - added support for new color list
% 2019-02-20: ABC - text changes with figure size

K = length(N);
color = get_color_list(K);

nfontsize2 = 0.25/3;


set(gca, 'XTick', [1:K]);
hold on
for i=1:K
  bar(i,N(i),'FaceColor', color{i});
end
title('ROI counts');
grid on
ylabel('count');
xlabel('ROI');

set(gca, 'FontUnits', 'normalized', 'FontSize', nfontsize2);
