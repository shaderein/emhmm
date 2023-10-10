function [data_new] = convert_face_coord(data, C, D, C_new, D_new)
% convert_face_coord - convert coordinates between face images.
%
%  [data_new] = convert_face_coord(data, C, D, C_new, D_new)
%
%  Convert the coordinate system from one face image to another face image 
%  of different size. This is helpful when using the standard analytic/holistic
%  models provided by the toolbox.
%  It uses the midpoint between the eyes and the distance from the 
%  eyes midline to the mouth to do the conversion.
% 
% INPUTS
%    data = the data to convert; data{i}{j} is the j-th trial of the i-th subject.
%       C = [x,y] coordinates of the midpoint between the eyes for the original image
%       D = distance from the eyes midline to the mouth for the original image
%   C_new = [x,y] coordinates of the midpoint between the eyes for the new image
%   D_new = distance from the eyes midline to the mouth for the new image
%
% OUTPUTS
%   data_new = the converted data for the new image size
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2017-05-25
% Antoni B. Chan, Janet H. Hsiao, Tim Chuk
% City University of Hong Kong, University of Hong Kong

N = length(data);
data_new = cell(1,N);

CC = C(:)';
CC_new = C_new(:)';

for i=1:N
  M = length(data{i});
  for j=1:M
    mydata = data{i}{j};
    tmp = bsxfun(@minus, mydata, CC)*(D_new/D);
    data_new{i}{j} = bsxfun(@plus, tmp, CC_new);
  end
end

