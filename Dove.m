clc, close, clear all;

%path/row- 170/078--J., South Africa
base='Z:\SpecialNeeds\Sakib\DoveData\P170\R078\';

dates = dir(base);
dates([1 2])=[];

equator= 10000000;
% ROI Coordinates
UL_x=433203;
UL_y=5730164-equator;
UR_x=433336;
UR_y=5730061-equator;
LR_x=433264;
LR_y=5729973-equator;
LL_x=433137;
LL_y=5730014-equator;

  
x_vec=[UL_x UR_x LR_x LL_x UL_x];
y_vec=[UL_y UR_y LR_y LL_y UL_y];
bands=  {'01' '02' '03' '04'}; 
band_name= {'Blue', 'Green', 'Red', 'NIR'};

%Image reading
Dove_image_file= strcat(base, dates.name, filesep, '35_q2.tif');
%Dove_image_file= strcat(base, filesep, '20181213_082751_35_1047_3B_AnalyticMS.tif');
[Dove_image_all_band, R_dove]=geotiffread(Dove_image_file);
% Dove_image_all_band=imread(Dove_image_file);
% R_dove =[0 -3; 3 0; 594918 7153524-equator ];

%Extracting all the values from Metadata file
MData_file=strcat(base, dates.name, filesep, '20181213_082751_35_1047_3B_AnalyticMS_metadata.xml');
[MData_values]= xml2struct_new_v(MData_file);

%Map coordinate to pixel coordinate
[Pixel_Row_unrounded_d, Pixel_Column_unrounded_d] = map2pix(R_dove, x_vec, y_vec);
Pixel_Row_d= round(Pixel_Row_unrounded_d);
Pixel_Column_d= round(Pixel_Column_unrounded_d);

% all the Angles
Dove_image_info.D_SEle_Angle=str2num(MData_values.ps_colon_EarthObservation.gml_colon_using...
    .eop_colon_EarthObservationEquipment.eop_colon_acquisitionParameters.ps_colon_Acquisition...
    .opt_colon_illuminationElevationAngle.Text);  
Dove_image_info.D_SAzi_Angle=str2num(MData_values.ps_colon_EarthObservation.gml_colon_using...
    .eop_colon_EarthObservationEquipment.eop_colon_acquisitionParameters.ps_colon_Acquisition...
    .opt_colon_illuminationAzimuthAngle.Text);  

% Loading Landsat 8 Data
load('L8_R.mat')
load('L8_ROIrad_b2.mat') 
load('L8_ROIrad_b3.mat')
load('L8_ROIrad_b4.mat')
load('L8_ROIrad_b5.mat')

 
% Blue Band-- Band 1 of Dove
for band=1:4
    Dove_image=Dove_image_all_band(:,:,band); %Blue band
    Dove_image=double(Dove_image);

    [Row_d, Column_d]= size(Dove_image);
    mask_d= poly2mask(Pixel_Column_d, Pixel_Row_d, Row_d, Column_d);
    Image_dove_masked=Dove_image.*mask_d;

    band_sf=0.01;
    D_ROIrad=Image_dove_masked.*(band_sf);
    %D_ROIrad=D_ROIrad./sind(Dove_image_info.D_SEle_Angle);
    D_ROIrad(D_ROIrad==0)=nan;
    %D_ROIrad_sampled=imresize(D_ROIrad, 3/30); % sampling

    % New corner Coordinates
    [nR,nC]=map2pix(R, R_dove.XWorldLimits(1,1), R_dove.YWorldLimits(1,2));
    nR=round(nR);
    nC=round(nC);
    [sR, sC]= size(D_ROIrad);
    
    D_ROIrad_final = NaN(size(L8_ROIrad_b2,1), size(L8_ROIrad_b2,2));%size is same for all L8 band
    D_ROIrad_final(nR:nR+sR-1,nC:nC+sC-1)=  D_ROIrad;

    if band==1
        D_ROIrand_band1=D_ROIrad_final;
    elseif band==2
        D_ROIrand_band2=D_ROIrad_final;
    elseif band==3
        D_ROIrand_band3=D_ROIrad_final;
    elseif band==4
        D_ROIrand_band4=D_ROIrad_final;
    end
end


%% Green band-- Band2 of Dove
Image_dove_b2=Dove_image(:,:,2);%Green Band
Image_dove_b2=double(Image_dove_b2);

[Row_g, Column_g]= size(Image_dove_b2);
mask_g= poly2mask(Pixel_Column2, Pixel_Row2, Row_g, Column_g);
Image_dove_b2_masked=Image_dove_b2.*mask_g;

TOArad_D_b2=Image_dove_b2_masked.*(0.01);
TOArad_D_b2=TOArad_D_b2./sind(37.41);
TOArad_D_b2(TOArad_D_b2==0)=nan;

TOArad_D_b2_sampled=imresize(TOArad_D_b2, 3/30); % sampling

load('L8_R.mat')
load('L8_ROIrad_b3.mat')

[nR,nC]=map2pix(R, 643431, 5020320);
nR=round(nR);
nC=round(nC);

TOArad_D_b2_final = NaN(size(L8_ROIrad_b3,1), size(L8_ROIrad_b3,2));
TOArad_D_b2_final(nR:nR+718-1,nC:nC+915-1)=TOArad_D_b2_sampled;


%Red band-- Band 3 of Dove
Image_dove_b3=Dove_image(:,:,3);%Red Band
Image_dove_b3=double(Image_dove_b3);

[Row_r, Column_r]= size(Image_dove_b3);
mask_r= poly2mask(Pixel_Column2, Pixel_Row2, Row_r, Column_r);
Image_dove_b3_masked=Image_dove_b3.*mask_r;

TOArad_D_b3=Image_dove_b3_masked.*(0.01);
TOArad_D_b3=TOArad_D_b3./sind(37.41);
TOArad_D_b3(TOArad_D_b3==0)=nan;

TOArad_D_b3_sampled=imresize(TOArad_D_b3, 3/30); % sampling

load('L8_R.mat')
load('L8_ROIrad_b4.mat')

[nR,nC]=map2pix(R, 643431, 5020320);
nR=round(nR);
nC=round(nC);

TOArad_D_b3_final = NaN(size(L8_ROIrad_b4,1), size(L8_ROIrad_b4,2));
TOArad_D_b3_final(nR:nR+718-1,nC:nC+915-1)=TOArad_D_b3_sampled;


%NIR band-- Band 4 of Dove
Image_dove_b4=Dove_image(:,:,4);%NIR Band
Image_dove_b4=double(Image_dove_b4);

[Row_nir, Column_nir]= size(Image_dove_b4);
mask_nir= poly2mask(Pixel_Column2, Pixel_Row2, Row_nir, Column_nir);
Image_dove_b4_masked=Image_dove_b4.*mask_nir;

TOArad_D_b4=Image_dove_b4_masked.*(0.01);
TOArad_D_b4=TOArad_D_b4./sind(37.41);
TOArad_D_b4(TOArad_D_b4==0)=nan;

TOArad_D_b4_sampled=imresize(TOArad_D_b4, 3/30); % sampling

load('L8_R.mat')
load('L8_ROIrad_b5.mat')

[nR,nC]=map2pix(R, 643431, 5020320);
nR=round(nR);
nC=round(nC);

TOArad_D_b4_final = NaN(size(L8_ROIrad_b5,1), size(L8_ROIrad_b5,2));
TOArad_D_b4_final(nR:nR+718-1,nC:nC+915-1)=TOArad_D_b4_sampled;

worldfile = getworldfilename(imagefile);
R = worldfileread(worldfile, 'geographic', size(RGB));

