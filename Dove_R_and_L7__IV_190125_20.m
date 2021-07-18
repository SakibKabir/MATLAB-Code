clearvars -except all_mean_DR all_SD_DR all_SNR_DR
bands=  {'1' '2' '3' '4' '5' '6' '7' '8'}; 
% %path/row- 
base='Z:\SpecialNeeds\Sakib\DoveData\April Data\Data\north_american_ag\imperial_valley_ag\P038\R038';
dates = dir(base); dates([1 2])=[];
date = '20190125';
Directory = dir(fullfile (base, date));
Directory([1 2]) = [];

%MTL File Parsing for Landsat 7
MTL=dir(fullfile(base, date, Directory.name,'*MTL.txt'));
[MTL_List_L7, vaule]= MTL_parser_L8(fullfile(MTL.folder, MTL.name));

% Dove File 
Dove_tif_name = '20190125_184150_20_1068_3B_AnalyticMS.tif';

%Extracting all the values from Metadata file
MData_file=fullfile(base, date, Directory.name, '20190125_184150_20_1068_3B_AnalyticMS_metadata.xml');
[MData_values]= xml2struct_new_v(MData_file);
[MData_values_old]= xml2struct(MData_file);

for ROI = 1: 4
  if ROI==1
    % ROI Coordinates
     UL_x= 692535;  UL_y= 3601637;
     UR_x= 692869;  UR_y= 3601564;
     LR_x= 692872;  LR_y= 3601424;
     LL_x= 692524;  LL_y= 3601522;

  elseif ROI==2
 UL_x= 698307;  UL_y= 3602737;
 UR_x= 698428;  UR_y= 3602840;
 LR_x= 698768;  LR_y= 3602450;
 LL_x= 698631;  LL_y= 3602372;

  elseif ROI==3
 UL_x= 707969;  UL_y= 3601296;
 UR_x= 708198;  UR_y= 3601294;
 LR_x= 708200;  LR_y= 3601163;
 LL_x= 707970;  LL_y= 3601160;

  elseif ROI==4
 UL_x= 689809;  UL_y= 3614601;
 UR_x= 689966;  UR_y= 3614736;
 LR_x= 690101;  LR_y= 3614460;
 LL_x= 689923;  LL_y= 3614339;

  end

 x_vec=[UL_x UR_x LR_x LL_x UL_x];
 y_vec=[UL_y UR_y LR_y LL_y UL_y];
  
for D_band=1:4 %dove band
    
L7_band = D_band;
L7_Folder_info=dir(fullfile(base, date, Directory.name, strcat(Directory.name,'_B', bands{L7_band},'.tif')));  
L7_image_file= fullfile(L7_Folder_info.folder, L7_Folder_info.name);
L7_image_info = geotiffinfo(L7_image_file);

[DN_L7, R_L7] = geotiffread(L7_image_file);
DN_L7 = double(DN_L7);
DN_L7(DN_L7 == 0)= nan;

% Mulitplicative and additive factors
rmb= MTL_List_L7.RADIOMETRIC_RESCALING.(strcat('RADIANCE_MULT_BAND_', bands{L7_band}));
rab= MTL_List_L7.RADIOMETRIC_RESCALING.(strcat('RADIANCE_ADD_BAND_', bands{L7_band}));
% L7 TOA Radiance
L7_TOArad = DN_L7*rmb+rab;

% Solar angles- scene center
Sun_Azimuth_L7 = MTL_List_L7.GROUP_IMAGE_ATTRIBUTES.SUN_AZIMUTH;
Sun_Elevation_L7 = MTL_List_L7.GROUP_IMAGE_ATTRIBUTES.SUN_ELEVATION;
L7_TOArad = L7_TOArad./sind(Sun_Elevation_L7); %cosine correction
MeanTOA_L7 = nanmean(nanmean(L7_TOArad));

%%% Landsat 7 TOA radiance of ROI
%Map coordinate to pixel coordinate
[Pixel_Row_unrounded_L7, Pixel_Column_unrounded_L7] = map2pix(R_L7, x_vec, y_vec);
Pixel_Row_L7= round(Pixel_Row_unrounded_L7);
Pixel_Column_L7= round(Pixel_Column_unrounded_L7);
[Row_L7, Column_L7]= size(DN_L7);

% ROI Mask 
mask_L7 = poly2mask(Pixel_Column_L7, Pixel_Row_L7,  Row_L7, Column_L7);
L7_TOArad_ROI = L7_TOArad.*mask_L7;
L7_TOArad_ROI(L7_TOArad_ROI == 0)= NaN;
MeanTOA_L7_ROI = nanmean(nanmean(L7_TOArad_ROI));

%%% Landsat 7 TOA reflectance of ROI
% Mulitplicative and additive factors
refmb= MTL_List_L7.RADIOMETRIC_RESCALING.(strcat('REFLECTANCE_MULT_BAND_', bands{L7_band}));
refab= MTL_List_L7.RADIOMETRIC_RESCALING.(strcat('REFLECTANCE_ADD_BAND_', bands{L7_band}));

L7_TOAref = DN_L7*refmb+refab;
L7_TOAref = L7_TOAref./sind(Sun_Elevation_L7); %cosine correction
MeanTOAref_L7 = nanmean(nanmean(L7_TOAref));

L7_TOAref_ROI = L7_TOAref.*mask_L7;
L7_TOAref_ROI(L7_TOAref_ROI == 0)= NaN;
MeanTOAref_L7_ROI = nanmean(nanmean(L7_TOAref_ROI));

%%% Dove 
Dove_Folder_info=dir(fullfile(base, date, Directory.name, Dove_tif_name));  
Dove_image_file= fullfile(Dove_Folder_info.folder, Dove_Folder_info.name);

Dove_image_all_band =imread(Dove_image_file);
DN_Dove = Dove_image_all_band(:,:, D_band);
DN_Dove = double(DN_Dove);
DN_Dove(DN_Dove == 0) = NaN;
band_sf=0.01;
Dove_TOArad = DN_Dove*band_sf;
MeanTOA_Dove = nanmean(nanmean(Dove_TOArad));

% Dove Cosine Correction
% all the Angles
Dove_image_info.D_SEle_Angle = str2double(MData_values.ps_colon_EarthObservation.gml_colon_using...
            .eop_colon_EarthObservationEquipment.eop_colon_acquisitionParameters.ps_colon_Acquisition...
            .opt_colon_illuminationElevationAngle.Text);  
Dove_image_info.D_SAzi_Angle = str2double(MData_values.ps_colon_EarthObservation.gml_colon_using...
            .eop_colon_EarthObservationEquipment.eop_colon_acquisitionParameters.ps_colon_Acquisition...
            .opt_colon_illuminationAzimuthAngle.Text);  
        
% Cosine Corrected Dove TOA Radiance inside the ROI        
Dove_TOArad = Dove_TOArad./sind(Dove_image_info.D_SEle_Angle);
MeanTOA_Dove = nanmean(nanmean(Dove_TOArad));

% Checking SNR and Standard Deviation image to select the ROI
window = 3;
Dove_mean_image = conv2(Dove_TOArad, ones(window)./(window*window),'same');
k= (window-1)/2; 
%k_L8= (window_L8-1)/2; 
Dove_SD_image = movingstd2(Dove_TOArad, k);
SNR_esti = Dove_mean_image./Dove_SD_image; 
% figure, imshow(Dove_SD_image, [])
% figure, imagesc(Dove_SD_image, [0 30]); colorbar
% figure, imagesc(Dove_SD_image); colorbar


% Dove Reference Matrix 
lat_TL = str2double(MData_values.ps_colon_EarthObservation.gml_colon_target...
            .ps_colon_Footprint.ps_colon_geographicLocation.ps_colon_topLeft.ps_colon_latitude.Text);   
lon_TL = str2double(MData_values.ps_colon_EarthObservation.gml_colon_target...
            .ps_colon_Footprint.ps_colon_geographicLocation.ps_colon_topLeft.ps_colon_longitude.Text);
[X, Y] = ll2utm(lat_TL, lon_TL);

RefMat_dove = [0 -3; 3 0; 686904 3620358];

%Map coordinate to pixel coordinate
[Pixel_Row_unrounded_d, Pixel_Column_unrounded_d] = map2pix(RefMat_dove, x_vec, y_vec);
Pixel_Row_d = round(Pixel_Row_unrounded_d);
Pixel_Column_d = round(Pixel_Column_unrounded_d);
[Row_d, Column_d] = size(DN_Dove);

mask_d= poly2mask(Pixel_Column_d, Pixel_Row_d, Row_d, Column_d);
Dove_TOArad_ROI = Dove_TOArad.*mask_d;
Dove_TOArad_ROI(Dove_TOArad_ROI == 0)= NaN;
MeanTOA_Dove_ROI = nanmean(nanmean(Dove_TOArad_ROI));

%%%% TOA Reflectance
%%% Reflectance conversion coefficient
    if D_band == 1
        id = 12; 
    elseif D_band == 2
        id = 14;
    elseif D_band == 3
        id = 16;
    elseif D_band == 4
        id = 18;
    end
    
Ref_coeff =str2double(MData_values_old.Children(10).Children(2).Children(id).Children(10).Children.Data);
Dove_TOAref = DN_Dove*Ref_coeff; 
MeanTOAref_Dove = nanmean(nanmean(Dove_TOAref));

% Cosine Corrected Dove TOA Reflectance inside the ROI 
Dove_TOAref = Dove_TOAref./sind(Dove_image_info.D_SEle_Angle);
Dove_TOAref_ROI = Dove_TOAref.*mask_d;
Dove_TOAref_ROI(Dove_TOAref_ROI == 0)= NaN;

%%% Radiance
temp_L7_rad = L7_TOArad_ROI;
temp_L7_rad=(temp_L7_rad(~isnan(temp_L7_rad)));

temp_D_rad = Dove_TOArad_ROI;
temp_D_rad=(temp_D_rad(~isnan(temp_D_rad)));

% Mean
mean_D_TOArad(D_band)= mean(temp_D_rad);
mean_L7_TOArad(D_band)= mean(temp_L7_rad); 

% Standard Deviation
Std_D_TOArad(D_band)= std(temp_D_rad);
Std_L7_TOArad(D_band)= std(temp_L7_rad);

%%% Reflectance
temp_D_ref = Dove_TOAref_ROI;
temp_D_ref=(temp_D_ref(~isnan(temp_D_ref)));

temp_L7_ref = L7_TOAref_ROI;
temp_L7_ref=(temp_L7_ref(~isnan(temp_L7_ref)));

% Mean
mean_D_TOAref(D_band)= mean(temp_D_ref);
mean_L7_TOAref(D_band)= mean(temp_L7_ref);

%Standard Deviation
Std_D_TOAref(D_band)= std(temp_D_ref);
Std_L7_TOAref(D_band)= std(temp_L7_ref);

end
 %%% Mean TOA radiance and reflectance
 %Mean- Radiance
 all_mean_IV_190125_20.D.ROI(ROI,:) = mean_D_TOArad;
 all_mean_IV_190125_20.L7.ROI(ROI,:) = mean_L7_TOArad;
 
 %Mean- Reflectance
 all_mean_ref_IV_190125_20.D.ROI(ROI,:) = mean_D_TOAref;
 all_mean_ref_IV_190125_20.L7.ROI(ROI,:) = mean_L7_TOAref;

 %%% Standard Deviation
 %Standard Deviation- Radiance
 all_SD_IV_190125_20.D.ROI(ROI,:)= Std_D_TOArad;
 all_SD_IV_190125_20.L7.ROI(ROI,:)= Std_L7_TOArad;

 %Standard Deviation- Reflectance
 all_SD_ref_IV_190125_20.D.ROI(ROI,:)= Std_D_TOAref;
 all_SD_ref_IV_190125_20.L7.ROI(ROI,:)= Std_L7_TOAref;

 %%% Signal to Noise Ratio
 %SNR From Radiance
 all_SNR_IV_190125_20.D.ROI(ROI,:) = mean_D_TOArad./Std_D_TOArad;
 all_SNR_IV_190125_20.L7.ROI(ROI,:) = mean_L7_TOArad./Std_L7_TOArad;

 %SNR From Reflectance
 all_SNR_ref_IV_190125_20.D.ROI(ROI,:) = mean_D_TOAref./Std_D_TOAref;
 all_SNR_ref_IV_190125_20.L7.ROI(ROI,:) = mean_L7_TOAref./Std_L7_TOAref; 
end

all_mean_DR.Radiance.a_IV_190125_20 = all_mean_IV_190125_20;
all_mean_DR.Reflectance.a_IV_190125_20 = all_mean_ref_IV_190125_20;

all_SD_DR.Radiance.a_IV_190125_20 = all_SD_IV_190125_20;
all_SD_DR.Reflectance.a_IV_190125_20= all_SD_ref_IV_190125_20;

all_SNR_DR.Radiance.a_IV_190125_20 = all_SNR_IV_190125_20;
all_SNR_DR.Reflectance.a_IV_190125_20 = all_SNR_ref_IV_190125_20;

save('all_mean_DR.mat','all_mean_DR')
save('all_SD_DR.mat',  'all_SD_DR')
save('all_SNR_DR.mat', 'all_SNR_DR')