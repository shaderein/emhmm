function len = get_median_length(data)
% get_median_length - get the median sequence length.
%
% len = get_median_length(data)
%
% INPUTS
%    data - cell array of sequences (one subject).
%         - or cell array of cell array of sequences (multiple subjects).
%         - (assumes that the sequence length is the first dimension)
% OUTPUT
%   len = median sequence length. 
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-08-12
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

% 2017-08-12: initial version
% 2018-11-23: v0.74 - add empty check

% get all lengths
L = get_lengths(data);

% calculate median
len = median(L);

% recursive function to get all lengths
function L = get_lengths(c)
if isempty(c)
  L = []; % do not include empty data
else
  if iscell(c{1})
    % it's a cell array of cells, so iterate on each cell
    L = [];
    for i=1:length(c)
      L = [L; get_lengths(c{i})];
    end
  else
    % it's data, so get the lengths
    L = cellfun(@(x) size(x,1), c);
    L = L(:);
  end
end
