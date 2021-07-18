clear all

for location = 1:19
% L8 and Dove Coincident file location and date extraction    
[base_L8, date_L8, base_D, date_D] = D_L8_Locations(location);

% For the southern hemisphere scene
switch location
    case{1,2,3,4,5,6,7,8}
        equator = 0;
    otherwise 
        equator = 10000000; 
end

% Dove Image File Path
DoveDataPath = fullfile(base_D, date_D, '*_AnalyticMS.tif');
DoveDirectory = dir(DoveDataPath);

% Dove XML File Path
DoveDataPathXML = fullfile(base_D, date_D, '*_metadata.xml');
DoveDirectoryXML = dir(DoveDataPathXML);
NoOf_Scene = size(DoveDirectory, 1); % total number of dove scene
%%
for scene_no = 1: NoOf_Scene 
    % Landsat 8 water filter from NIR band, give percentage as decimal
    percentage = 0.05;
    
    L8_w_filter = waterEli_L8_Ref(base_L8, date_L8, percentage);
    % figure, imagesc(L8_w_filter); colorbar
    
    D_band_filter = 4;  percentage_D = 0.08;
    Dove_w_filter = waterEli_D_Ref(base_D, D_band_filter, scene_no, percentage_D); 
    % figure, imagesc(Dove_w_filter); colorbar
    %%
% for band = 1:4
    %% Calculating Landsat 8 TOA Reflectance
    % Blue Band
    [L8_TOAref_B, R_L8] = TOAref_cal_L8(base_L8, date_L8, 2);
    % figure, imagesc(L8_TOAref_B); colorbar
    % figure, histogram(L8_TOAref_B)
    % Green Band
    L8_TOAref_G = TOAref_cal_L8(base_L8, date_L8, 3);
    % Red Band
    L8_TOAref_R = TOAref_cal_L8(base_L8, date_L8, 4);
    % NIR Band
    L8_TOAref_NIR = TOAref_cal_L8(base_L8, date_L8, 5);
    
    %%%% Applying water filter to Landsat 8 Scene
    L8_TOAref_B = L8_TOAref_B.*L8_w_filter; clear L8_w_filter
    % figure, imagesc(L8_TOAref_Wfil); colorbar
    
    %%%% Calculating Dove TOA Reflectance
    %D_band = band - 1; % equator for southern and northern hemisphere scene
    [Dove_TOAref_B, R_Dove] = TOAref_cal_D(base_D, scene_no, 1, equator);
    % figure, imagesc(Dove_TOAref_B); colorbar
    Dove_TOAref_G = TOAref_cal_D(base_D, scene_no, 2, equator);
    Dove_TOAref_R = TOAref_cal_D(base_D, scene_no, 3, equator);
    Dove_TOAref_NIR = TOAref_cal_D(base_D, scene_no, 4, equator);
    
    %%%% Applying water filter to Dove Scene
    Dove_TOAref_B = Dove_TOAref_B.*Dove_w_filter; clear Dove_w_filter
    % figure, imagesc(Dove_TOAref); colorbar
      
    %% Size of the ROI 
    window = 15;  
    k = (window-1)/2; 
    
    %%% L8 Mean image
    L8_mean_image = conv2(L8_TOAref_B, ones(window)./(window*window),'same'); clear window
    % figure, imagesc(L8_mean_image); colorbar
    % figure, histogram(L8_mean_image)
    
    %%% L8 Standard Deviation image
    L8_SD_image = movingstd2(L8_TOAref_B, k); clear k
    % figure, imagesc(L8_SD_image); colorbar
    % figure, histogram(L8_SD_image)

    %%% Coefficient of Variation 
    L8_CV = L8_SD_image./L8_mean_image; %clear L8_mean_image L8_SD_image
    % figure, imagesc(L8_CV); colorbar
    % figure, histogram(L8_CV)
    
    %%% ROI Selection 
    L8_CV_ROI_15 = L8_CV; clear L8_CV
    CV_low = 0; CV_high = 0.01; % As a start: don't take it as a threshold
    L8_CV_ROI_15(L8_CV_ROI_15 < CV_low | L8_CV_ROI_15 > CV_high) = NaN; 
    % figure, imagesc(L8_CV_ROI_15); colorbar
    % figure, histogram(L8_CV_ROI_15)
    clear CV_low CV_high
    
    L8_CV_ROI_15_temp = L8_CV_ROI_15;
    L8_CV_ROI_15_temp(~isnan(L8_CV_ROI_15_temp)) = 1;
    % figure, imagesc(L8_CV_ROI_15_temp); colorbar
    
    L8_TOAref_temp = L8_TOAref_B;
    L8_TOAref_temp = L8_CV_ROI_15_temp.*L8_TOAref_temp;
    % figure, imagesc(L8_TOAref_temp); colorbar
    clear L8_CV_ROI_15 L8_CV_ROI_15_temp 
    
    %%
    % Creating an OLI size NaN matrix
    L8_nan_mat = NaN(size(L8_TOAref_B, 1), size(L8_TOAref_B, 2));
    % Creating Mean, SD and CV NaN matrix to store L8 values
    L8_mean_image_ROI5 = L8_nan_mat;
    L8_SD_image_ROI5 = L8_nan_mat; clear L8_nan_mat
   
    % Creating a Dove size NaN matrix
    D_nan_mat = nan(size(Dove_TOAref_B, 1), size(Dove_TOAref_B, 2));
    % Creating Mean, SD and CV NaN matrix to store Dove values
    D_mean_image_ROI5 = D_nan_mat; 
    D_SD_image_ROI5 = D_nan_mat; clear D_nan_mat

    % Extracting all the values from Metadata file
    MData_file = fullfile(base_D, date_D, DoveDirectoryXML(scene_no).name);
    [MData_values]= xml2struct_new_v(MData_file);
    
    % Lat Lon Top Left
    lat_TL = str2num(MData_values.ps_colon_EarthObservation.gml_colon_target...
             .ps_colon_Footprint.ps_colon_geographicLocation.ps_colon_topLeft.ps_colon_latitude.Text);   
    lon_TL = str2num(MData_values.ps_colon_EarthObservation.gml_colon_target...
            .ps_colon_Footprint.ps_colon_geographicLocation.ps_colon_topLeft.ps_colon_longitude.Text);
    
    % Lat Lon Top Right    
    lat_TR = str2num(MData_values.ps_colon_EarthObservation.gml_colon_target...
            .ps_colon_Footprint.ps_colon_geographicLocation.ps_colon_topRight.ps_colon_latitude.Text);   
    lon_TR = str2num(MData_values.ps_colon_EarthObservation.gml_colon_target...
            .ps_colon_Footprint.ps_colon_geographicLocation.ps_colon_topRight.ps_colon_longitude.Text);
        
    % Lat Lon Bottom Right    
    lat_BR = str2num(MData_values.ps_colon_EarthObservation.gml_colon_target...
            .ps_colon_Footprint.ps_colon_geographicLocation.ps_colon_bottomRight.ps_colon_latitude.Text);   
    lon_BR = str2num(MData_values.ps_colon_EarthObservation.gml_colon_target...
            .ps_colon_Footprint.ps_colon_geographicLocation.ps_colon_bottomRight.ps_colon_longitude.Text);
        
    % Lat Lon Bottom Left    
    lat_BL = str2num(MData_values.ps_colon_EarthObservation.gml_colon_target...
            .ps_colon_Footprint.ps_colon_geographicLocation.ps_colon_bottomLeft.ps_colon_latitude.Text);   
    lon_BL = str2num(MData_values.ps_colon_EarthObservation.gml_colon_target...
            .ps_colon_Footprint.ps_colon_geographicLocation.ps_colon_bottomLeft.ps_colon_longitude.Text);
    clear MData_values MData_file
    
    lat = [lat_TL lat_TR lat_BR lat_BL]; clear lat_TL lat_TR lat_BR lat_BL
    lon = [lon_TL lon_TR lon_BR lon_BL]; clear lon_TL lon_TR lon_BR lon_BL
    % Lat-lon to UTM conversion
    [X, Y] = ll2utm(lat, lon); clear lat lon
    % map to pixel conversion L8
    [L8_pixel_row, L8_pixel_col] = map2pix(R_L8, X, Y-equator); clear X Y
    L8_pixel_row = round(L8_pixel_row); 
    L8_pixel_col = round(L8_pixel_col);
    [Row_L8, Column_L8]= size(L8_TOAref_B);

    % ROI Mask 
    mask_L8 = poly2mask(L8_pixel_col, L8_pixel_row,  Row_L8, Column_L8); clear L8_pixel_col L8_pixel_row
    L8_TOAref_ROI = L8_TOAref_temp.*mask_L8; clear L8_TOAref_temp mask_L8
    L8_TOAref_ROI(L8_TOAref_ROI == 0)= NaN;
    % figure, imagesc(L8_TOAref_ROI); colorbar
    clear Row_L8 Column_L8
    % finding the rows and columns of the pixel where CV is in the range
    [rows, cols] = find(~isnan(L8_TOAref_ROI)); %clear L8_TOAref_ROI
    
    %%
for r = 1: size(rows, 1)
    %% Rows and columns are in pixel cordinates
     temp_rows = rows(r);
     temp_cols = cols(r);
%     temp_rows = 5701;
%     temp_cols = 5545;

    % Extracting the matrix and taking the mean and sd for 5 by 5
    % pixel from L8_TOAref matrix
    L8_TOAref_ROI_5B = L8_TOAref_B(temp_rows-2:temp_rows+2, temp_cols-2:temp_cols+2);
    L8_meanTOAref_ROI_5B(r) = nanmean(L8_TOAref_ROI_5B(:)); 
    L8_sdTOAref_ROI_5B(r) = nanstd(L8_TOAref_ROI_5B(:)); 
    
    L8_TOAref_ROI_5G = L8_TOAref_G(temp_rows-2:temp_rows+2, temp_cols-2:temp_cols+2);
    L8_meanTOAref_ROI_5G(r) = nanmean(L8_TOAref_ROI_5G(:)); 
    L8_sdTOAref_ROI_5G(r) = nanstd(L8_TOAref_ROI_5G(:)); clear L8_TOAref_ROI_5G
    
    L8_TOAref_ROI_5R = L8_TOAref_R(temp_rows-2:temp_rows+2, temp_cols-2:temp_cols+2);
    L8_meanTOAref_ROI_5R(r) = nanmean(L8_TOAref_ROI_5R(:)); 
    L8_sdTOAref_ROI_5R(r) = nanstd(L8_TOAref_ROI_5R(:)); clear L8_TOAref_ROI_5R
    
    L8_TOAref_ROI_5NIR = L8_TOAref_NIR(temp_rows-2:temp_rows+2, temp_cols-2:temp_cols+2);
    L8_meanTOAref_ROI_5NIR(r) = nanmean(L8_TOAref_ROI_5NIR(:)); 
    L8_sdTOAref_ROI_5NIR(r) = nanstd(L8_TOAref_ROI_5NIR(:)); clear L8_TOAref_ROI_5NIR
    
    %%%% Pasting mean, sd and cv at exact OLI location to visually inspect for sanity check
    L8_mean_image_ROI5(temp_rows, temp_cols) = nanmean(L8_TOAref_ROI_5B(:));
    % L8_SD_image_ROI5(temp_rows, temp_cols) = nanstd(L8_TOAref_ROI_5(:));
    % L8_CV_image_ROI5(temp_rows, temp_cols) = nanstd(L8_TOAref_ROI_5(:))/nanmean(L8_TOAref_ROI_5(:));
    % figure, imagesc(L8_mean_image_ROI5); colorbar
    
    % Getting the map cordinates from OLI pixel cordinates for Dove 
    temp_row_vec = [temp_rows-2, temp_rows-2, temp_rows+3, temp_rows+3]; clear temp_rows
    temp_col_vec = [temp_cols-2, temp_cols+3, temp_cols+3, temp_cols-2]; clear temp_cols
    [x_vec, y_vec] = pix2map(R_L8, temp_row_vec, temp_col_vec); clear temp_row_vec temp_col_vec
    
    %%%% Dove Overlaping Regions TOA reflectance calculation
    % Map(from L8) to Pixel cordinate
    [Pixel_Row_unrounded_D, Pixel_Column_unrounded_D] = map2pix(R_Dove, x_vec, y_vec);
    Pixel_Row_D = round(Pixel_Row_unrounded_D); 
    Pixel_Column_D = round(Pixel_Column_unrounded_D);
    clear x_vec y_vec Pixel_Row_unrounded_D Pixel_Column_unrounded_D
    %%
    if (Pixel_Row_D > 1 & Pixel_Row_D <=  size(Dove_TOAref_B, 1)) & (Pixel_Column_D > 1 &...
            Pixel_Column_D <=  size(Dove_TOAref_B, 2))
       
       % Dove Mean and Sd Calculation
       Dove_mean_TOA_ROI_B = nanmean(nanmean(Dove_TOAref_B(min(Pixel_Row_D):max(Pixel_Row_D),...
                             min(Pixel_Column_D):max(Pixel_Column_D))));
       Dove_sd_TOA_ROI_B = nanstd(nanstd(Dove_TOAref_B(min(Pixel_Row_D):max(Pixel_Row_D),...
                           min(Pixel_Column_D):max(Pixel_Column_D)))); 
       
       % Green Band
       Dove_mean_TOA_ROI_G = nanmean(nanmean(Dove_TOAref_G(min(Pixel_Row_D):max(Pixel_Row_D),...
                             min(Pixel_Column_D):max(Pixel_Column_D))));
       Dove_sd_TOA_ROI_G = nanstd(nanstd(Dove_TOAref_G(min(Pixel_Row_D):max(Pixel_Row_D),...
                             min(Pixel_Column_D):max(Pixel_Column_D)))); 
       % Red Band
       Dove_mean_TOA_ROI_R = nanmean(nanmean(Dove_TOAref_R(min(Pixel_Row_D):max(Pixel_Row_D),...
                             min(Pixel_Column_D):max(Pixel_Column_D))));
       Dove_sd_TOA_ROI_R = nanstd(nanstd(Dove_TOAref_R(min(Pixel_Row_D):max(Pixel_Row_D),...
                             min(Pixel_Column_D):max(Pixel_Column_D)))); 
       
       % NIR Band
       Dove_mean_TOA_ROI_NIR = nanmean(nanmean(Dove_TOAref_NIR(min(Pixel_Row_D):max(Pixel_Row_D),...
                             min(Pixel_Column_D):max(Pixel_Column_D))));
       Dove_sd_TOA_ROI_NIR = nanstd(nanstd(Dove_TOAref_NIR(min(Pixel_Row_D):max(Pixel_Row_D),...
                             min(Pixel_Column_D):max(Pixel_Column_D)))); 
                       
       % Storing all the values
       Dove_meanTOAref_ROI_5B(r) = Dove_mean_TOA_ROI_B;
       Dove_sdTOAref_ROI_5B(r) = Dove_sd_TOA_ROI_B;
       
       Dove_meanTOAref_ROI_5G(r) = Dove_mean_TOA_ROI_G; clear Dove_mean_TOA_ROI_G
       Dove_sdTOAref_ROI_5G(r) = Dove_sd_TOA_ROI_G; clear Dove_sd_TOA_ROI_G
       Dove_meanTOAref_ROI_5R(r) = Dove_mean_TOA_ROI_R; clear Dove_mean_TOA_ROI_R
       Dove_sdTOAref_ROI_5R(r) = Dove_sd_TOA_ROI_R; clear Dove_sd_TOA_ROI_R
       Dove_meanTOAref_ROI_5NIR(r) = Dove_mean_TOA_ROI_NIR; clear Dove_mean_TOA_ROI_NIR
       Dove_sdTOAref_ROI_5NIR(r) = Dove_sd_TOA_ROI_NIR; clear Dove_sd_TOA_ROI_NIR
       
       
       % Inserting mean, sd  at exact dove pixel location for sanity check
       % Middle row and column location 
       D_mean_location_row = (min(Pixel_Row_D) + max(Pixel_Row_D))/2; 
       D_mean_location_col = (min(Pixel_Column_D) + max(Pixel_Column_D))/2;
       clear Pixel_Row_D Pixel_Column_D
       
       % creation of matrix
       D_mean_image_ROI5(D_mean_location_row, D_mean_location_col) = Dove_mean_TOA_ROI_B; clear Dove_mean_TOA_ROI_B
       % figure, imagesc(D_mean_image_ROI5), colorbar
       [r1, c1] = find(D_mean_image_ROI5 == max(max(D_mean_image_ROI5)));
       D_SD_image_ROI5(D_mean_location_row, D_mean_location_col) = Dove_sd_TOA_ROI_B; clear Dove_sd_TOA_ROI_B
       clear D_mean_location_row D_mean_location_col Dove_mean_TOA_ROI Dove_sd_TOA_ROI
       
    else
       Dove_meanTOAref_ROI_5B(r) = nan;
       Dove_sdTOAref_ROI_5B(r) = nan;
       
       Dove_meanTOAref_ROI_5G(r) = nan;
       Dove_sdTOAref_ROI_5G(r) = nan;
       Dove_meanTOAref_ROI_5R(r) = nan;
       Dove_sdTOAref_ROI_5R(r) = nan;
       Dove_meanTOAref_ROI_5NIR(r) = nan;
       Dove_sdTOAref_ROI_5NIR(r) = nan;
       
    end

end; clear r Pixel_Row_D Pixel_Column_D

   if size(rows, 1) == 0
          % To keep the structure of the matrix consistent
%           L8_meanTOAref_ROI(:,band) = {[]}; L8_sdTOAref_ROI(:,band) = {[]};
%           D_meanTOAref_ROI(:,band) = {[]}; D_sdTOAref_ROI(:,band) = {[]};
     D_meanTOAref_ROI = []; D_sdTOAref_ROI = []; L8_meanTOAref_ROI = [];
     L8_sdTOAref_ROI = [];
   else

     D_meanTOAref_ROI = [Dove_meanTOAref_ROI_5B; Dove_meanTOAref_ROI_5G; Dove_meanTOAref_ROI_5R; Dove_meanTOAref_ROI_5NIR]';
     D_sdTOAref_ROI = [Dove_sdTOAref_ROI_5B; Dove_sdTOAref_ROI_5G; Dove_sdTOAref_ROI_5R; Dove_sdTOAref_ROI_5NIR]';
     clear Dove_meanTOAref_ROI_5B Dove_meanTOAref_ROI_5G Dove_meanTOAref_ROI_5R Dove_meanTOAref_ROI_5NIR
     clear Dove_sdTOAref_ROI_5B Dove_sdTOAref_ROI_5G Dove_sdTOAref_ROI_5R Dove_sdTOAref_ROI_5NIR

     L8_meanTOAref_ROI= [L8_meanTOAref_ROI_5B; L8_meanTOAref_ROI_5G; L8_meanTOAref_ROI_5R; L8_meanTOAref_ROI_5NIR]';
     L8_sdTOAref_ROI= [L8_sdTOAref_ROI_5B; L8_sdTOAref_ROI_5G; L8_sdTOAref_ROI_5R; L8_sdTOAref_ROI_5NIR]';
     clear L8_meanTOAref_ROI_5B L8_meanTOAref_ROI_5G L8_meanTOAref_ROI_5R L8_meanTOAref_ROI_5NIR
     clear L8_sdTOAref_ROI_5B L8_sdTOAref_ROI_5G L8_sdTOAref_ROI_5R L8_sdTOAref_ROI_5NIR
 
   end 
   % Storing all the data
   L8_mTOAref_allB(scene_no, 1) = {L8_meanTOAref_ROI}; clear L8_meanTOAref_ROI
   L8_sdTOAref_allB(scene_no, 1) = {L8_sdTOAref_ROI}; clear L8_sdTOAref_ROI
   D_mTOAref_allB(scene_no, 1) = {D_meanTOAref_ROI}; clear D_meanTOAref_ROI
   D_sdTOAref_allB(scene_no, 1) = {D_sdTOAref_ROI}; clear D_sdTOAref_ROI
   
end
clear date dates dates_L8 equator base_D base_L8 
     CCTOAref.D_mean(location, 1) = {D_mTOAref_allB}; clear D_mTOAref_allB
     CCTOAref.L8_mean(location, 1) = {L8_mTOAref_allB}; clear L8_mTOAref_allB
     CCTOAref.D_sd(location, 1) = {D_sdTOAref_allB}; clear D_sdTOAref_allB
     CCTOAref.L8_sd(location, 1) = {L8_sdTOAref_allB}; clear L8_sdTOAref_allB
end
clear Dove_TOAref_B Dove_TOAref_G Dove_TOAref_R Dove_TOAref_NIR
clear L8_TOAref_B L8_TOAref_G L8_TOAref_R L8_TOAref_NIR
clearvars -except CCTOAref
save('CCTOAref2.mat')