function S = print_cell_table(T, width)
% print_cell_table - convert a cell array of values into a nice table for display
%
%   S = print_cell_table(T, width)
%
% INPUTS
%       T - cell array
%   width - width for each column
%       S - output string
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2019-07-13
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

% 2020-03-12: v0.77 - initial version


S = '';

for i=1:size(T,1)
  for j=1:size(T,2)
    if isinteger(T{i,j})
      tmp = sprintf('%d', T{i,j});
      
    elseif isnumeric(T{i,j})
      tmp = sprintf('%0.4g', T{i,j});
      
    else
      tmp = T{i,j};
    end
    S = [S pad(tmp, width, 'left')];
    S = [S ' '];
  end
  if i < size(T,1)
    S = [S sprintf('\n')];
  end
end


function a = pad(b, width, lr)

diff = width - length(b);
if (diff>0)
  switch(lr)
    case 'left'
      a = [repmat(' ', 1, diff), b];
    case 'right'
      a = [b, repmat(' ', 1, diff)];
    otherwise
      a = b;     
  end
else
  a = b;
end
