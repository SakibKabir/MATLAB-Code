function [L8_TOAref, R_L8] = TOAref_cal_L8(base_L8, date_L8, band)
    %% Landsat 8
    bands =  {'1' '2' '3' '4' '5' };
    
    %%% MTL File Parsing 
    MTL = dir(fullfile(base_L8, date_L8, 'LC1', '*MTL.txt'));
    [MTL_List, vaule]= MTL_parser_L8(fullfile(MTL.folder, MTL.name));
    
    % Image file location
    Band_No = dir(fullfile(base_L8, date_L8, 'LC1', strcat('*B', bands{band},'.TIF')));
    L8_image_file = fullfile(base_L8, date_L8, 'LC1', Band_No.name);

    % BQA File location and Reading the BQA Image, creating mask
    BQA_File_dir = dir(fullfile(base_L8, date_L8, 'LC1', '*_BQA.TIF'));
    L8_BQA_File = fullfile(base_L8, date_L8, 'LC1', BQA_File_dir.name);
    [BQA_data, R_BQA] = geotiffread(L8_BQA_File);
    BQA_mask =(BQA_data == 2720);

    % Reading the OLI image
    [DN_L8, R_L8] = geotiffread(L8_image_file);
    DN_L8 = double(DN_L8);
    DN_L8(DN_L8 == 0)= nan;
    L8_image_info = geotiffinfo(L8_image_file);
    % figure, imagesc(DN_L8); colorbar

    % Mulitplicative and additive factors
    refmb= MTL_List.RADIOMETRIC_RESCALING.(strcat('REFLECTANCE_MULT_BAND_', bands{band}));
    refab= MTL_List.RADIOMETRIC_RESCALING.(strcat('REFLECTANCE_ADD_BAND_', bands{band}));

    % L8 TOA Radiance
    L8_TOAref_raw = DN_L8*refmb+refab;
    % figure, histogram(L8_TOArad_raw); 

    % L8 TOA Radiance masked
    L8_TOAref = L8_TOAref_raw.*BQA_mask;
    L8_TOAref(L8_TOAref == 0)= nan;
    % figure, imagesc(L8_TOArad); colorbar
    % figure, histogram(L8_TOArad);

    % L8 TOA Radiance cosine correction: Scene Center Angle
    % Sun_Elevation = MTL_List.GROUP_IMAGE_ATTRIBUTES.SUN_ELEVATION;
    % L8_TOArad = L8_TOArad./sind(Sun_Elevation); %cosine correction
    % figure, imagesc(L8_TOArad); colorbar
    % figure, histogram(L8_TOArad);

    % Angles for cosine correction: Pixel wise
    Angle_folder= dir(fullfile(base_L8, date_L8, 'LC1', strcat('*solar_B05.img')));
    L8_Angle_file= fullfile(Angle_folder.folder, Angle_folder.name);
    SolarAngleInfo = multibandread(L8_Angle_file,[size(DN_L8),2],'int16',0,'bsq','ieee-le');
    %solar_azimuth=(SolarAngleInfo(:,:,1))/100;
    solar_zenith = (SolarAngleInfo(:,:,2))/100;


    %L8 solar zenith pixel by pixel
    L8_mat_logical = ~isnan(L8_TOAref);
    solar_zenith_L8ROImat = solar_zenith.*L8_mat_logical;
    solar_zenith_L8ROImat(solar_zenith_L8ROImat==0) = nan;

    L8_TOAref = L8_TOAref./cosd(solar_zenith_L8ROImat);  
    % figure, imagesc(L8_TOArad); colorbar
    % figure, histogram(L8_TOArad);
end 