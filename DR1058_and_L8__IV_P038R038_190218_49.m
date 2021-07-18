%clearvars -except all_mean_DR_1058 all_SD_DR_1058 all_SNR_DR_1058
bands=  {'1' '2' '3' '4' '5' '6' '7' '8'}; 
% %path/row- 
base='Z:\SpecialNeeds\Sakib\DoveData\April Data\Dove-R Data\1058\20190218\P038\R038';
dates = dir(base); dates([1 2])=[];
date = '20190218';
Directory = dir(fullfile (base, date));
Directory([1 2]) = [];

%MTL File Parsing for Landsat 7
MTL=dir(fullfile(base, date, Directory.name,'*MTL.txt'));
[MTL_List_L8, vaule]= MTL_parser_L8(fullfile(MTL.folder, MTL.name));

% Dove File 
Dove_tif_name = '20190218_180350_49_1058_3B_AnalyticMS.tif';

%Extracting all the values from Metadata file
MData_file=fullfile(base, date, Directory.name, '20190218_180350_49_1058_3B_AnalyticMS_metadata.xml');
[MData_values]= xml2struct_new_v(MData_file);
[MData_values_old]= xml2struct(MData_file);

for ROI = 1: 4
  if ROI==1
    % ROI Coordinates
 UL_x= 705551;  UL_y= 3565028;
 UR_x= 706025;  UR_y= 3565062;
 LR_x= 706106;  LR_y= 3564338;
 LL_x= 705555;  LL_y= 3564417;

  elseif ROI==2
 UL_x= 691295;  UL_y= 3568488;
 UR_x= 691551; UR_y= 3568491;
 LR_x= 691555;  LR_y= 3568335;
 LL_x= 691305;  LL_y= 3568333;

  elseif ROI==3
%  UL_x= 693390;  UL_y= 3566874;
%  UR_x= 693496; UR_y= 3566815;
%  LR_x= 693343;  LR_y= 3566711;
%  LL_x= 693257;  LL_y= 3566781;

 UL_x= 688639;  UL_y= 3554499;
 UR_x= 688670;  UR_y= 3554476;
 LR_x= 688655;  LR_y= 3554452;
 LL_x= 688629;  LL_y= 3554473;
 
  elseif ROI==4
%  UL_x= 689766;  UL_y= 3558703;
%  UR_x= 689959;  UR_y= 3558509;
%  LR_x= 689767;  LR_y= 3558340;
%  LL_x= 689651;  LL_y= 3558441;

 UL_x= 694479;  UL_y= 3566292;
 UR_x= 694566;  UR_y= 3566328;
 LR_x= 694649;  LR_y= 3566191;
 LL_x= 694564;  LL_y= 3566148;
  end

 x_vec=[UL_x UR_x LR_x LL_x UL_x];
 y_vec=[UL_y UR_y LR_y LL_y UL_y];
  
for D_band=1:4 %dove band
    
L8_band = D_band+1;
L8_Folder_info=dir(fullfile(base, date, Directory.name, strcat(Directory.name,'_B', bands{L8_band},'.tif')));  
L8_image_file= fullfile(L8_Folder_info.folder, L8_Folder_info.name);
L8_image_info = geotiffinfo(L8_image_file);

[DN_L8, R_L8] = geotiffread(L8_image_file);
DN_L8 = double(DN_L8);
DN_L8(DN_L8 == 0)= nan;

% Mulitplicative and additive factors
rmb= MTL_List_L8.RADIOMETRIC_RESCALING.(strcat('RADIANCE_MULT_BAND_', bands{L8_band}));
rab= MTL_List_L8.RADIOMETRIC_RESCALING.(strcat('RADIANCE_ADD_BAND_', bands{L8_band}));
% L8 TOA Radiance
L8_TOArad = DN_L8*rmb+rab;

% Solar angles- scene center
Sun_Azimuth_L8 = MTL_List_L8.GROUP_IMAGE_ATTRIBUTES.SUN_AZIMUTH;
Sun_Elevation_L8 = MTL_List_L8.GROUP_IMAGE_ATTRIBUTES.SUN_ELEVATION;
L8_TOArad = L8_TOArad./sind(Sun_Elevation_L8); %cosine correction
MeanTOA_L8 = nanmean(nanmean(L8_TOArad));

%%% Landsat 7 TOA radiance of ROI
%Map coordinate to pixel coordinate
[Pixel_Row_unrounded_L8, Pixel_Column_unrounded_L8] = map2pix(R_L8, x_vec, y_vec);
Pixel_Row_L8= round(Pixel_Row_unrounded_L8);
Pixel_Column_L8= round(Pixel_Column_unrounded_L8);
[Row_L8, Column_L8]= size(DN_L8);

% ROI Mask 
mask_L8 = poly2mask(Pixel_Column_L8, Pixel_Row_L8,  Row_L8, Column_L8);
L8_TOArad_ROI = L8_TOArad.*mask_L8;
L8_TOArad_ROI(L8_TOArad_ROI == 0)= NaN;
MeanTOA_L8_ROI = nanmean(nanmean(L8_TOArad_ROI));

%%% Landsat 7 TOA reflectance of ROI
% Mulitplicative and additive factors
refmb= MTL_List_L8.RADIOMETRIC_RESCALING.(strcat('REFLECTANCE_MULT_BAND_', bands{L8_band}));
refab= MTL_List_L8.RADIOMETRIC_RESCALING.(strcat('REFLECTANCE_ADD_BAND_', bands{L8_band}));

L8_TOAref = DN_L8*refmb+refab;
L8_TOAref = L8_TOAref./sind(Sun_Elevation_L8); %cosine correction
MeanTOAref_L8 = nanmean(nanmean(L8_TOAref));

L8_TOAref_ROI = L8_TOAref.*mask_L8;
L8_TOAref_ROI(L8_TOAref_ROI == 0)= NaN;
MeanTOAref_L8_ROI = nanmean(nanmean(L8_TOAref_ROI));

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
% figure, imagesc(Dove_SD_image, [0 10]); colorbar
% figure, imagesc(Dove_SD_image); colorbar

% Dove Reference Matrix 
lat_TL = str2double(MData_values.ps_colon_EarthObservation.gml_colon_target...
            .ps_colon_Footprint.ps_colon_geographicLocation.ps_colon_topLeft.ps_colon_latitude.Text);   
lon_TL = str2double(MData_values.ps_colon_EarthObservation.gml_colon_target...
            .ps_colon_Footprint.ps_colon_geographicLocation.ps_colon_topLeft.ps_colon_longitude.Text);
[X, Y] = ll2utm(lat_TL, lon_TL);

RefMat_dove = [0 -3; 3 0; 687342 3569691];

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
%Dove_TOAref = Dove_TOAref./sind(Dove_image_info.D_SEle_Angle);
Dove_TOAref_ROI = Dove_TOAref.*mask_d;
Dove_TOAref_ROI(Dove_TOAref_ROI == 0)= NaN;

%%% Radiance
temp_L8_rad = L8_TOArad_ROI;
temp_L8_rad=(temp_L8_rad(~isnan(temp_L8_rad)));

temp_D_rad = Dove_TOArad_ROI;
temp_D_rad=(temp_D_rad(~isnan(temp_D_rad)));

% Mean
mean_D_TOArad(D_band)= mean(temp_D_rad);
mean_L8_TOArad(D_band)= mean(temp_L8_rad); 

% Standard Deviation
Std_D_TOArad(D_band)= std(temp_D_rad);
Std_L8_TOArad(D_band)= std(temp_L8_rad);

%%% Reflectance
temp_D_ref = Dove_TOAref_ROI;
temp_D_ref=(temp_D_ref(~isnan(temp_D_ref)));

temp_L8_ref = L8_TOAref_ROI;
temp_L8_ref=(temp_L8_ref(~isnan(temp_L8_ref)));

% Mean
mean_D_TOAref(D_band)= mean(temp_D_ref);
mean_L8_TOAref(D_band)= mean(temp_L8_ref);

%Standard Deviation
Std_D_TOAref(D_band)= std(temp_D_ref);
Std_L8_TOAref(D_band)= std(temp_L8_ref);

end
 %%% Mean TOA radiance and reflectance
 %Mean- Radiance
 all_mean_IV_P038R038_190218_49.D.ROI(ROI,:) = mean_D_TOArad;
 all_mean_IV_P038R038_190218_49.L8.ROI(ROI,:) = mean_L8_TOArad;
 
 %Mean- Reflectance
 all_mean_ref_IV_P038R038_190218_49.D.ROI(ROI,:) = mean_D_TOAref;
 all_mean_ref_IV_P038R038_190218_49.L8.ROI(ROI,:) = mean_L8_TOAref;

 %%% Standard Deviation
 %Standard Deviation- Radiance
 all_SD_IV_P038R038_190218_49.D.ROI(ROI,:)= Std_D_TOArad;
 all_SD_IV_P038R038_190218_49.L8.ROI(ROI,:)= Std_L8_TOArad;

 %Standard Deviation- Reflectance
 all_SD_ref_IV_P038R038_190218_49.D.ROI(ROI,:)= Std_D_TOAref;
 all_SD_ref_IV_P038R038_190218_49.L8.ROI(ROI,:)= Std_L8_TOAref;

 %%% Signal to Noise Ratio
 %SNR From Radiance
 all_SNR_IV_P038R038_190218_49.D.ROI(ROI,:) = mean_D_TOArad./Std_D_TOArad;
 all_SNR_IV_P038R038_190218_49.L8.ROI(ROI,:) = mean_L8_TOArad./Std_L8_TOArad;

 %SNR From Reflectance
 all_SNR_ref_IV_P038R038_190218_49.D.ROI(ROI,:) = mean_D_TOAref./Std_D_TOAref;
 all_SNR_ref_IV_P038R038_190218_49.L8.ROI(ROI,:) = mean_L8_TOAref./Std_L8_TOAref; 
end

% all_mean_DR_1058.Radiance.a_IV_P038R038_190218_49 = all_mean_IV_P038R038_190218_49;
% all_mean_DR_1058.Reflectance.a_IV_P038R038_190218_49 = all_mean_ref_IV_P038R038_190218_49;
% 
% all_SD_DR_1058.Radiance.a_IV_P038R038_190218_49 = all_SD_IV_P038R038_190218_49;
% all_SD_DR_1058.Reflectance.a_IV_P038R038_190218_49= all_SD_ref_IV_P038R038_190218_49;
% 
% all_SNR_DR_1058.Radiance.a_IV_P038R038_190218_49 = all_SNR_IV_P038R038_190218_49;
% all_SNR_DR_1058.Reflectance.a_IV_P038R038_190218_49 = all_SNR_ref_IV_P038R038_190218_49;
% 
% save('all_mean_DR_1058.mat','all_mean_DR_1058')
% save('all_SD_DR_1058.mat',  'all_SD_DR_1058')
% save('all_SNR_DR_1058.mat', 'all_SNR_DR_1058')