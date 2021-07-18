clc, close, clear all;

%path/row- 170/078--J., South Africa
base='Z:\SpecialNeeds\Sakib\DoveData\P170\R078';
dates = dir(base);
dates([1 2])=[];

equator= 10000000;
% ROI Coordinates
UL_x=594920;
UL_y=7148999-equator;
UR_x=619285;
UR_y=7153522-equator;
LR_x=622482;
LR_y=7137102-equator;
LL_x=598082;
LL_y=7132591-equator;

x_vec=[UL_x UR_x LR_x LL_x UL_x];
y_vec=[UL_y UR_y LR_y LL_y UL_y];

bands=  {'1' '2' '3' '4' '5' '6' '7'}; 
band_str= {'CA','B',  'G', 'R', 'NIR', 'SWIR1', 'SWIR2'};

%MTL File Parsing 
MTL=dir(fullfile(base, dates.name,'*MTL.txt'));
[MTL_List, vaule]= MTL_parser_L8(fullfile(MTL.folder, MTL.name));

for band=2:5
    %Image file location
    Folder_info=dir(fullfile(base, dates.name, strcat('*B', bands{band},'.tif')));  
    L8_image_file= fullfile(Folder_info.folder, Folder_info.name);

    %Reading the image
    [DN, R] = geotiffread(L8_image_file);
    DN=double(DN);
    L8_info = geotiffinfo(L8_image_file);

    %Map coordinate to pixel coordinate
    [Pixel_Row_unrounded, Pixel_Column_unrounded] = map2pix(R, x_vec, y_vec);
    Pixel_Row= round(Pixel_Row_unrounded);
    Pixel_Column= round(Pixel_Column_unrounded);
    [Row, Column]= size(DN);

    % ROI Mask 
    mask= poly2mask(Pixel_Column, Pixel_Row,  Row, Column);
    DN=DN.*mask;

    % Mulitplicative and additive factors
    rmb= MTL_List.RADIOMETRIC_RESCALING.(strcat('RADIANCE_MULT_BAND_', bands{band}));
    rab= MTL_List.RADIOMETRIC_RESCALING.(strcat('RADIANCE_ADD_BAND_', bands{band}));
    L8_ROIrad=DN*rmb+rab;

    % Solar angles
    Sun_Azimuth= MTL_List.GROUP_IMAGE_ATTRIBUTES.SUN_AZIMUTH;
    Sun_Elevation=MTL_List.GROUP_IMAGE_ATTRIBUTES.SUN_ELEVATION;

    %L8_ROIrad_b2=L8_ROIrad_b2./sind(Sun_Elevation); %cosine correction
    L8_ROIrad(L8_ROIrad<0)=nan;

    %Storing radiances to different varibles
    if band==2
        L8_ROIrad_b2=L8_ROIrad;
        clear L8_ROIrad
    elseif band==3
        L8_ROIrad_b3=L8_ROIrad;
        clear L8_ROIrad
    elseif band==4
        L8_ROIrad_b4=L8_ROIrad;
        clear L8_ROIrad
    elseif band==5
        L8_ROIrad_b5=L8_ROIrad;
        clear L8_ROIrad
    end
end


%%
%GREEN band
band=3;
L8_image_file= strcat(base, filesep, 'LC08_L1TP_026029_20180923_20180929_01_T1_B3.tif');
[DN, R] = geotiffread(L8_image_file);
DN=double(DN);
L8_info = geotiffinfo(L8_image_file);

%Map coordinate to pixel coordinate
[Pixel_Row_unrounded, Pixel_Column_unrounded] = map2pix(R, x_vec, y_vec);
Pixel_Row= round(Pixel_Row_unrounded);
Pixel_Column= round(Pixel_Column_unrounded);
[Row, Column]= size(DN);

% ROI Mask 
mask= poly2mask(Pixel_Column, Pixel_Row,  Row, Column);
DN=DN.*mask;

% Mulitplicative and additive factors
rmb= MTL_List.RADIOMETRIC_RESCALING.(strcat('RADIANCE_MULT_BAND_', bands{band}));
rab= MTL_List.RADIOMETRIC_RESCALING.(strcat('RADIANCE_ADD_BAND_', bands{band}));
L8_ROIrad_b2=DN*rmb-rab;

L8_ROIrad_b3=DN*0.011769-58.84354;
L8_ROIrad_b3=L8_ROIrad_b3./sind(42.7);
L8_ROIrad_b3(L8_ROIrad_b3<0)=nan;


 %RED band
L8_image_file= strcat(base, filesep, 'LC08_L1TP_026029_20180923_20180929_01_T1_B4.tif');
[DN, R] = geotiffread(L8_image_file);
DN=double(DN);
L8_info = geotiffinfo(L8_image_file);

%Map coordinate to pixel coordinate
[Pixel_Row_unrounded, Pixel_Column_unrounded] = map2pix(R, x_vec, y_vec);
Pixel_Row= round(Pixel_Row_unrounded);
Pixel_Column= round(Pixel_Column_unrounded);
[Row, Column]= size(DN);

% ROI Mask 
mask= poly2mask(Pixel_Column, Pixel_Row,  Row, Column);
DN=DN.*mask;
L8_ROIrad_b4=DN*0.009924-49.62019;
L8_ROIrad_b4=L8_ROIrad_b4./sind(42.7);
L8_ROIrad_b4(L8_ROIrad_b4<0)=nan;

% NIR Band
L8_image_file= strcat(base, filesep, 'LC08_L1TP_026029_20180923_20180929_01_T1_B5.tif');
[DN, R] = geotiffread(L8_image_file);
DN=double(DN);
L8_info = geotiffinfo(L8_image_file);

%Map coordinate to pixel coordinate
[Pixel_Row_unrounded, Pixel_Column_unrounded] = map2pix(R, x_vec, y_vec);
Pixel_Row= round(Pixel_Row_unrounded);
Pixel_Column= round(Pixel_Column_unrounded);
[Row, Column]= size(DN);

% ROI Mask 
mask= poly2mask(Pixel_Column, Pixel_Row,  Row, Column);
DN=DN.*mask;
L8_ROIrad_b5=DN*0.006073-30.36508;
L8_ROIrad_b5=L8_ROIrad_b5./sind(42.7);
L8_ROIrad_b5(L8_ROIrad_b5<0)=nan;
