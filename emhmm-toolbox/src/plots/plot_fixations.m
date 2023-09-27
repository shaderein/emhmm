function plot_fixations(data, img, LL, opt)
% plot_fixations - plot fixations on an image
%
%  plot_fixations(data, img, LL, opt)
%
% INPUT: 
%   data  = data (same format as used by vbhmm_learn)
%   img   = image for background
%   LL    = the log-likelihood (for display)
%           [] - to ignore it
%   opt   = 's' - show sequences [default] 
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-01-13
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong

% VERSIONS
% 2017-01-17: ABC - add duration as marker size.
% 2017-05-25: ABC - added options to not show sequence lines
% 2017-06-21: ABC - use our own "showim" instead of "imshow"

if nargin<2
  img = '';
end
if nargin<3
  LL = [];
end
if nargin<4
  opt = 's';
end

nfontsize2 = 0.25/3;

if ~isempty(data)
  D = size(data{1}, 2);
  if D > 2
    myscale = get_duration_rescale(data);
  end
end

if ~isempty(img)
  showim(img);
end
axis ij
hold on
for i=1:length(data)
  if (D > 2)
    % plot duration as marker size
    szs1 = myscale(data{i}(1,3));
    szs2 = myscale(data{i}(2:end,3));
  else
    % fixed marker size
    szs1 = 6;
    szs2 = 6;
  end
  % plot x,y location, and duration as marker size
  if any(opt=='s')
    plot(data{i}(:,1), data{i}(:,2), 'b-');
  end
  scatter(data{i}(1,1), data{i}(1,2), szs1, 'co'); % first fixation
  scatter(data{i}(2:end,1), data{i}(2:end,2), szs2, 'bo'); % others
end
hold off
if ~isempty(LL)
  title(sprintf('fixations (LL=%g)', LL));
  set(gca, 'FontUnits', 'normalized', 'FontSize', nfontsize2);
end

