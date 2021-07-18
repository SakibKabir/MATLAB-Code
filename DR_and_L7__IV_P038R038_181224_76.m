clearvars -except all_mean_DR_105e all_SD_DR_105e all_SNR_DR_105e
bands=  {'1' '2' '3' '4' '5' '6' '7' '8'}; 
% %path/row- 
base='Z:\SpecialNeeds\Sakib\DoveData\April Data\Dove-R Data\105e\20181224\P038\R038';
dates = dir(base); dates([1 2])=[];
date = '20181224';
Directory = dir(fullfile (base, date));
Directory([1 2]) = [];

%MTL File Parsing for Landsat 7
MTL=dir(fullfile(base, date, Directory.name,'*MTL.txt'));
[MTL_List_L7, vaule]= MTL_parser_L8(fullfile(MTL.folder, MTL.name));

% Dove File 
Dove_tif_name = '20181224_180028_76_105e_3B_AnalyticMS.tif';

%Extracting all the values from Metadata file
MData_file=fullfile(base, date, Directory.name, '20181224_180028_76_105e_3B_AnalyticMS_metadata.xml');
[MData_values]= xml2struct_new_v(MData_file);
[MData_values_old]= xml2struct(MData_file);

for ROI = 1:4
  if ROI==1
    % ROI Coordinates
 UL_x= 694004;  UL_y= 3556700;
 UR_x= 694141;  UR_y= 3556696;
 LR_x= 694148;  LR_y= 3556522;
 LL_x= 693990;  LL_y= 3556442;

  elseif ROI==2
 UL_x= 694082;  UL_y= 3560202;
 UR_x= 694494;  UR_y= 3560062;
 LR_x= 694452;  LR_y= 3559868;

 LL_x1= 694344;  LL_y1= 3559890;
 LL_x2= 694085;  LL_y2= 3560145;

LL_x = [LL_x1 LL_x2];
LL_y=  [LL_y1 LL_y2];

  elseif ROI==3
 UL_x= 712033;  UL_y= 3553075;
 UR_x= 712859;  UR_y= 3552981;
 LR_x= 713093;  LR_y= 3552322;
 LL_x= 712071;  LL_y= 3552402;

  elseif ROI==4
 UL_x= 697155;  UL_y= 3555901;
 UR_x= 697214;  UR_y= 3555908;

 LR_x1= 697252;  LR_y1= 3555864;
 LR_x2= 697355;  LR_y2= 3555887;
 LR_x3= 697382;  LR_y3= 3555832;
 
 LR_x = [LR_x1 LR_x2 LR_x3];
 LR_y=  [LR_y1 LR_y2 LR_y3];

 LL_x1= 697314;  LL_y1= 3555781;
 LL_x2= 697266;  LL_y2= 3555800;
 LL_x3= 697192;  LL_y3= 3555796;

LL_x = [LL_x1 LL_x2 LL_x3];
LL_y=  [LL_y1 LL_y2 LL_y3];

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
L7_TOAref_ROI(L7_TOAref_ROI < 0)= NaN;
MeanTOAref_L7_ROI = nanmean(nanmean(L7_TOAref_ROI));

%%% Dove 
Dove_Folder_info=dir(fullfile(base, date, Directory.name, Dove_tif_name));  
Dove_image_file= fullfile(Dove_Folder_info.folder, Dove_Folder_info.name);
%Dove_image_info = geotiffinfo(Dove_image_file);

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
% figure, imagesc(Dove_SD_image, [0 15]); colorbar
% figure, imagesc(Dove_SD_image); colorbar


% Dove Reference Matrix 
lat_TL = str2double(MData_values.ps_colon_EarthObservation.gml_colon_target...
            .ps_colon_Footprint.ps_colon_geographicLocation.ps_colon_topLeft.ps_colon_latitude.Text);   
lon_TL = str2double(MData_values.ps_colon_EarthObservation.gml_colon_target...
            .ps_colon_Footprint.ps_colon_geographicLocation.ps_colon_topLeft.ps_colon_longitude.Text);
[X, Y] = ll2utm(lat_TL, lon_TL);

RefMat_dove = [0 -3; 3 0; 690714 3562041];

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
 all_mean_IV_P038R038_181224_76.D.ROI(ROI,:) = mean_D_TOArad;
 all_mean_IV_P038R038_181224_76.L7.ROI(ROI,:) = mean_L7_TOArad;
 
 %Mean- Reflectance
 all_mean_ref_IV_P038R038_181224_76.D.ROI(ROI,:) = mean_D_TOAref;
 all_mean_ref_IV_P038R038_181224_76.L7.ROI(ROI,:) = mean_L7_TOAref;

 %%% Standard Deviation
 %Standard Deviation- Radiance
 all_SD_IV_P038R038_181224_76.D.ROI(ROI,:)= Std_D_TOArad;
 all_SD_IV_P038R038_181224_76.L7.ROI(ROI,:)= Std_L7_TOArad;

 %Standard Deviation- Reflectance
 all_SD_ref_IV_P038R038_181224_76.D.ROI(ROI,:)= Std_D_TOAref;
 all_SD_ref_IV_P038R038_181224_76.L7.ROI(ROI,:)= Std_L7_TOAref;

 %%% Signal to Noise Ratio
 %SNR From Radiance
 all_SNR_IV_P038R038_181224_76.D.ROI(ROI,:) = mean_D_TOArad./Std_D_TOArad;
 all_SNR_IV_P038R038_181224_76.L7.ROI(ROI,:) = mean_L7_TOArad./Std_L7_TOArad;

 %SNR From Reflectance
 all_SNR_ref_IV_P038R038_181224_76.D.ROI(ROI,:) = mean_D_TOAref./Std_D_TOAref;
 all_SNR_ref_IV_P038R038_181224_76.L7.ROI(ROI,:) = mean_L7_TOAref./Std_L7_TOAref; 
end

all_mean_DR_105e.Radiance.a_IV_P038R038_181224_76 = all_mean_IV_P038R038_181224_76;
all_mean_DR_105e.Reflectance.a_IV_P038R038_181224_76 = all_mean_ref_IV_P038R038_181224_76;

all_SD_DR_105e.Radiance.a_IV_P038R038_181224_76 = all_SD_IV_P038R038_181224_76;
all_SD_DR_105e.Reflectance.a_IV_P038R038_181224_76= all_SD_ref_IV_P038R038_181224_76;

all_SNR_DR_105e.Radiance.a_IV_P038R038_181224_76 = all_SNR_IV_P038R038_181224_76;
all_SNR_DR_105e.Reflectance.a_IV_P038R038_181224_76 = all_SNR_ref_IV_P038R038_181224_76;

save('all_mean_DR_105e.mat','all_mean_DR_105e')
save('all_SD_DR_105e.mat',  'all_SD_DR_105e')
save('all_SNR_DR_105e.mat', 'all_SNR_DR_105e')