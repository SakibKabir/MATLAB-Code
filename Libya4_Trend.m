clc, close, clear all;

base='Z:\ImageDrive\OLI.TIRS\L8\P181\R040';
dates = dir(base); dates([1 2])=[];
bands=  {'1' '2' '3' '4' '5' '6' '7'}; 

% CEOS ROI
UL_x=723825; UL_y=3171375;
UR_x=743355; UR_y=3171825;
LR_x=724245; LR_y=3149325;
LL_x=743805; LL_y=3149685;

ULlat= 3171375 ; ULlong= 723825;   %mapping coordinates
URlat= 3171825 ; URlong= 743355 ;
LRlat= 3149325; LRlong = 724245 ;
LLlat= 3149685; LLlong= 743805 ; 

x_vec=[UL_x UR_x LR_x LL_x UL_x];
y_vec=[UL_y UR_y LR_y LL_y UL_y];

%%
for date = 1:10%size(dates,1)
 for band = 1:7
%%
    Folder_info = dir(fullfile(base, dates(date).name, 'LC1', strcat('*B', bands{band},'.tif')));
    L8_image_file= fullfile(Folder_info.folder, Folder_info.name);
    
    %Reading the image
    [DN, R_L8] = geotiffread(L8_image_file);
    DN = double(DN);
    L8_info = geotiffinfo(L8_image_file);

    % MTL Parser Function to extract all the data from the Meta Data File
    MTL=dir(fullfile(base, dates(date).name,'LC1','*MTL.txt'));
    [MTL_List, vaule]= MTL_parser_L8(fullfile(MTL.folder, MTL.name));
    SE = MTL_List.GROUP_IMAGE_ATTRIBUTES.(strcat('SUN_ELEVATION'));
    
    % Day of Year and Decimal Year
    acquisition_date = MTL_List.PRODUCT_METADATA.DATE_ACQUIRED;      %Date extraction from MetaData file
    [doy,fraction] = date2doy(datenum(acquisition_date));            %date2doy function to convert date to get Day of the Year
    DoY(date) = doy;                                                  %day of year
    DateVec = datevec(acquisition_date);                              %converting acquisition date to date vector
    DeciYear(date)=DateVec(1,1)+DateVec(1,2)./12+DateVec(1,3)./365;  %decimal year
    
    rmb= MTL_List.RADIOMETRIC_RESCALING.REFLECTANCE_MULT_BAND_1;
    rab= MTL_List.RADIOMETRIC_RESCALING.REFLECTANCE_ADD_BAND_1;
    
    %Map coordinate to pixel coordinate
    [Pixel_Row_unrounded_L8, Pixel_Column_unrounded_L8] = map2pix(R_L8, x_vec, y_vec);
    Pixel_Row_L8= round(Pixel_Row_unrounded_L8);
    Pixel_Column_L8= round(Pixel_Column_unrounded_L8);
    [Row_L8, Column_L8]= size(DN);

    % ROI Mask 
    mask_L8= poly2mask(Pixel_Column_L8, Pixel_Row_L8,  Row_L8, Column_L8);
    DN=DN.*mask_L8;
    DN(DN==0)= NaN;
    
    % BQA File location and Reading the BQA Image, creating mask
    BQA_File_dir = dir(fullfile(base, dates(date).name, 'LC1', '*_BQA.TIF'));
    L8_BQA_File = strcat(base, filesep, dates(date).name, filesep, 'LC1', filesep, BQA_File_dir.name);
    [BQA_data, R] = geotiffread(L8_BQA_File);
    BQA_mask =(BQA_data == 2720);
    
    L8_TOAref = DN*rmb+rab;
    L8_TOAref(L8_TOAref == 0) = NaN;
    L8_TOAref = L8_TOAref.*BQA_mask;
    
    % Angles for pixelwise cosine correction
    Angle_folder= dir(fullfile(base, dates(date).name,'LC1', strcat('*solar_B05.img')));
    L8_Angle_file= fullfile(Angle_folder.folder, Angle_folder.name);
    SolarAngleInfo = multibandread(L8_Angle_file,[size(DN),2],'int16',0,'bsq','ieee-le');
    % solar_azimuth=(SolarAngleInfo(:,:,1))/100;
    solar_zenith =(SolarAngleInfo(:,:,2))/100;
    
    % L8 solar zenith pixel by pixel
    L8_mat_logical = ~isnan(L8_TOAref);
    solar_zenith_L8ROImat = solar_zenith.*L8_mat_logical;
    solar_zenith_L8ROImat(solar_zenith_L8ROImat==0) = nan;
    L8_TOAref_CC = L8_TOAref./cosd(solar_zenith_L8ROImat); 
    
    %%% Reflectance
    temp_D_ref = L8_TOAref_CC;
    temp_D_ref=(temp_D_ref(~isnan(temp_D_ref)));
    
    Libya4_CEOSref(date, band) = nanmean(L8_TOAref_CC(:));
    Libya4_SD(date, band) = nanstd(nanstd(L8_TOAref_CC));
 end
end
