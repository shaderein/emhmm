function showim(img)
% draw an image
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-06-21
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong

% 2018-05-16: v0.72 - load images from filename
% 2019-02-22: v0.75 - support grayscale mode

global EMHMM_COLORPLOTS

if ischar(img)
  imgname = img;
  img = imread(imgname);
end

% process the image
if EMHMM_COLORPLOTS == 0
  if ndims(img) == 3
    % v0.75 - convert color to grayscale
    dimg = double(img);
    gimg = 0.2989*dimg(:,:,1) + 0.5870*dimg(:,:,2) + 0.1140*dimg(:,:,3);
    img = uint8(gimg);
  end
end
  
if ndims(img) == 2
  img = cat(3, img, img, img);
end

image(img);
axis image;
axis off;