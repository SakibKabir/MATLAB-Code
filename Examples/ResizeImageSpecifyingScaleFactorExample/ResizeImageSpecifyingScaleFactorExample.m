%% Resize Image Specifying Scale Factor
% 
%%
% Read image into the workspace.

% Copyright 2015 The MathWorks, Inc.

I = imread('rice.png');
%%
% Resize the image, specifying scale factor and using default interpolation
% method and antialiasing.
J = imresize(I, [128 130]);
%%
% Display the original and the resized image.
figure
imshow(I)
title('Original Image')
figure
imshow(J)
title('Resized Image')
