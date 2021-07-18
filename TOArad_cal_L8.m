function [L8_TOArad, R_L8] = TOArad_cal_L8(base_L8, date, band)
    %% Landsat 8
    base_L8 = base_L8; date = date; band = band;
    dates_L8 = dir(base_L8); dates_L8([1 2])=[];
    bands =  {'1' '2' '3' '4' '5' }; 
    
    %%% MTL File Parsing 
    MTL = dir(fullfile(base_L8, dates_L8(date).name,'LC1','*MTL.txt'));
    [MTL_List, vaule]= MTL_parser_L8(fullfile(MTL.folder, MTL.name));

    % Image file location
    Band_No = dir(fullfile(base_L8, dates_L8(date).name, 'LC1', strcat('*B', bands{band},'.tif')));
    L8_image_file = strcat(base_L8, filesep, dates_L8(date).name, filesep, 'LC1', filesep, Band_No.name);

    % BQA File location and Reading the BQA Image, creating mask
    BQA_File_dir = dir(fullfile(base_L8, dates_L8(date).name, 'LC1', '*_BQA.TIF'));
    L8_BQA_File = strcat(base_L8, filesep, dates_L8(date).name, filesep, 'LC1', filesep, BQA_File_dir.name);
    [BQA_data, R] = geotiffread(L8_BQA_File);
    BQA_mask =(BQA_data == 2720);

    % Reading the OLI image
    [DN_L8, R_L8] = geotiffread(L8_image_file);
    DN_L8 = double(DN_L8);
    DN_L8(DN_L8 == 0)= nan;
    L8_image_info = geotiffinfo(L8_image_file);
    % figure, imagesc(DN_L8); colorbar

    % Mulitplicative and additive factors
    rmb = MTL_List.RADIOMETRIC_RESCALING.(strcat('RADIANCE_MULT_BAND_', bands{band}));
    rab = MTL_List.RADIOMETRIC_RESCALING.(strcat('RADIANCE_ADD_BAND_', bands{band}));

    % L8 TOA Radiance
    L8_TOArad_raw = DN_L8*rmb+rab;
    % figure, histogram(L8_TOArad_raw); 

    % L8 TOA Radiance masked
    L8_TOArad = L8_TOArad_raw.*BQA_mask;
    L8_TOArad(L8_TOArad == 0)= nan;
    % figure, imagesc(L8_TOArad); colorbar
    % figure, histogram(L8_TOArad);

    % L8 TOA Radiance cosine correction: Scene Center Angle
    % Sun_Elevation = MTL_List.GROUP_IMAGE_ATTRIBUTES.SUN_ELEVATION;
    % L8_TOArad = L8_TOArad./sind(Sun_Elevation); %cosine correction
    % figure, imagesc(L8_TOArad); colorbar
    % figure, histogram(L8_TOArad);

    % Angles for cosine correction: Pixel wise
    Angle_folder= dir(fullfile(base_L8, dates_L8(date).name, 'LC1', strcat('*solar_B05.img')));
    L8_Angle_file= fullfile(Angle_folder.folder, Angle_folder.name);
    SolarAngleInfo = multibandread(L8_Angle_file,[size(DN_L8),2],'int16',0,'bsq','ieee-le');
    %solar_azimuth=(SolarAngleInfo(:,:,1))/100;
    solar_zenith = (SolarAngleInfo(:,:,2))/100;


    %L8 solar zenith pixel by pixel
    L8_mat_logical = ~isnan(L8_TOArad);
    solar_zenith_L8ROImat = solar_zenith.*L8_mat_logical;
    solar_zenith_L8ROImat(solar_zenith_L8ROImat==0) = nan;

    L8_TOArad = L8_TOArad./cosd(solar_zenith_L8ROImat);  
    % figure, imagesc(L8_TOArad); colorbar
    % figure, histogram(L8_TOArad);
end 