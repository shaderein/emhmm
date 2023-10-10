function [outliers] = find_degenerate_outliers(all_LL, all_h3ms, mode)
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2019-07-13
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

% 2020-03-12: v0.77 - initial version

if nargin<3
  mode = 'diff';
end

switch(mode)
  case '3std'
    % method 1: use standard deviation
    mn = mean(all_LL);
    st = std(all_LL);
    outliers = find(all_LL > mn+2*st);

  case 'diff'
    % method 2: look at differences
    [all_LL_sort, ii] = sort(all_LL, 'descend');
    df = diff(all_LL_sort);

    outliers = [];
    i=1;
    % remove items that are > 10% better than the next best.
    while( abs(df(i) / all_LL_sort(i+1)) > 0.10 )
      outliers(end+1) = ii(i);
      i = i + 1;
    end
end


