function plot_face_coord(title1, faceimg1, C1, D1, title2, faceimg2, C2, D2)
% plot_face_coord - plot face coordinates for conversion
%
%  plot_face_coord(title1, faceimg1, C1, D1, title2, faceimg2, C2, D2)
%
% 
% INPUTS
%   title1, faceimg1, C1, D1 - 1st set of title, image, C and D 
%   title2, faceimg2, C2, D2 - 2nd set of title, image, C and D
%
% SEE ALSO
%    convert_face_coord
%
% ---
% Eye-Movement analysis with HMMs (emhmm-toolbox)
% Copyright (c) 2018-05-15
% Antoni B. Chan, Janet H. Hsiao
% City University of Hong Kong, University of Hong Kong

% VERSIONS
% 2018-05-16: v0.72 - initial version

figure
subplot(1,2,1)
showim(faceimg1);
hold on
plot(C1(1), C1(2), 'ro');
plot(C1(1)*[0.25, 1.75], C1(2)*[1 1], 'r-');
plot(C1(1), C1(2)+D1, 'bx');
plot(C1(1)*[1 1], C1(2)+[0, D1], 'b-');
hold off
legend({'midpoint between eyes', 'midline between eyes', ...
  'mouth level', 'distance from eyes midpoint to mouth level'},...
  'Location', 'SouthOutside');
title(title1)

subplot(1,2,2)
showim(faceimg2);
hold on
plot(C2(1), C2(2), 'ro');
plot(C2(1)*[0.25, 1.75], C2(2)*[1 1], 'r-');
plot(C2(1), C2(2)+D2, 'bx');
plot(C2(1)*[1 1], C2(2)+[0, D2], 'b-');
hold off
title(title2);
