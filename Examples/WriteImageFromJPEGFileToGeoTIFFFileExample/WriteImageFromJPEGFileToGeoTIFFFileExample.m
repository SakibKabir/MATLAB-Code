%% Write Image from JPEG File to GeoTIFF File
% 
%%
% Read JPEG image from file. 

% Copyright 2015 The MathWorks, Inc.
base='Z:\SpecialNeeds\Sakib\DoveData\P026\R029\LC08_L1TP_026029_20180923_20180929_01_T1';
Dove_image_file= strcat(base, filesep, '20180923_155545_1047_3B_AnalyticMS.tif');

basename = 'boston_ovr';
imagefile = [basename '.jpg'];
RGB = imread(Dove_image_file);
%%
% Derive world file name from image file name, read the world file, and
% construct a spatial referencing object. 
worldfile = getworldfilename('20180923_155545_1047_3B_AnalyticMS.tif');

R = worldfileread(worldfile, 'geographic', size(RGB));
%%
% Write image data and referencing data to GeoTIFF file.
filename = [basename '.tif'];
geotiffwrite(filename, RGB, R)
%%
% Construct an empty map axes and display the map.
figure
usamap(RGB, R)
geoshow(filename)
