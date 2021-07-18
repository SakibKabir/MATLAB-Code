clc, close, clear all;

base= 'Z:\ImageDrive\OLI.TIRS\L8\P181\R040';
%base= 'Z:\SpecialNeeds\Sakib\DoveData\P026\029';

dates = dir(base);
dates([1 2])=[];

% ROI Coordinates
UL_x=739080;
UL_y=3162691;
UR_x=770212;
UR_y=3162576;
LR_x=770498;
LR_y=3126914;
LL_x=731970;
LL_y=3128004;

x_vec=[UL_x UR_x LR_x LL_x UL_x];
y_vec=[UL_y UR_y LR_y LL_y UL_y];

bands=  {'1' '2' '3' '4' '5' '6' '7'}; 
band_str= {'CA','B',  'G', 'R', 'NIR', 'SWIR1', 'SWIR2'};

for band = 2 %:7
 for date = 52 %1 :size(dates,1)
%    for date = 1:3
   if exist (fullfile(base, dates(date).name, 'LC1'),'dir')
       
        %MTL File Parsing 
        MTL=dir(fullfile(base, dates(date).name,'LC1','*MTL.txt'));
        [MTL_List, vaule]= MTL_parser_L8(fullfile(MTL.folder, MTL.name));
       
        %Image file location
        Band_No(band)=dir(fullfile(base, dates(date).name, 'LC1', strcat('*B', bands{band},'.tif')));
        L8_image_file= strcat(base, filesep, dates(date).name, filesep, 'LC1', filesep, Band_No(band).name);
        
        % Mulitplicative and additive factors
        rmb(date)= MTL_List.RADIOMETRIC_RESCALING.(strcat('RADIANCE_MULT_BAND_', bands{band}));
        rab(date)= MTL_List.RADIOMETRIC_RESCALING.(strcat('RADIANCE_ADD_BAND_', bands{band}));
       
        % Mulitplicative and additive factors- reflectance for checking
        refmb(date)= MTL_List.RADIOMETRIC_RESCALING.(strcat('REFLECTANCE_MULT_BAND_', bands{band}));
        refab(date)= MTL_List.RADIOMETRIC_RESCALING.(strcat('REFLECTANCE_ADD_BAND_', bands{band}));
        
        %Reading Image
        [DN, R] = geotiffread(L8_image_file);
        DN=double(DN);
        L8_info=geotiffinfo(L8_image_file);
       
        %Map coordinate to pixel coordinate
        [Pixel_Row_unrounded, Pixel_Column_unrounded] = map2pix(R, x_vec, y_vec);
        Pixel_Row= round(Pixel_Row_unrounded);
        Pixel_Column= round(Pixel_Column_unrounded);
        [Row, Column]= size(DN);
        
        % ROI Mask 
        mask= poly2mask(Pixel_Column, Pixel_Row,  Row, Column);
        DN=DN.*mask;
        
        %TOA Radiance calculation
        TOArad_L8=DN*rmb(date)+rab(date); 
        TOArad_L8(TOArad_L8<0)=NaN;
       
        % Solar angles
        Sun_Azimuth(date)= MTL_List.GROUP_IMAGE_ATTRIBUTES.SUN_AZIMUTH;
        Sun_Elevation(date)=MTL_List.GROUP_IMAGE_ATTRIBUTES.SUN_ELEVATION;
        
        TOArad_L8=TOArad_L8./sind(Sun_Elevation(date));%TOA Radiance after cosine correction
        TOArad_L8(TOArad_L8<0)=NaN;
        
%         %%%checking
%         %TOA Reflectance- for checking only
%         TOAref_L8=DN*refmb(date)+refab(date);
%         TOAref_L8=TOAref_L8./sind(Sun_Elevation(date));
%         TOAref_L8( TOAref_L8<0)=NaN;
%   
        %MeanTOArad_L8(date) =mean(TOArad_L8); %Mean TOA Radiance
        
        % Day of Year and Decimal Year
        acquisition_date =MTL_List.PRODUCT_METADATA.DATE_ACQUIRED;
        [doy,fraction] = date2doy(datenum(acquisition_date));
        DoY(date)=doy; %day of year
        DateVec= datevec(acquisition_date);
        DeciYear(date)=DateVec(1,1)+DateVec(1,2)./12+DateVec(1,3)./365; %decimal year
   
   end 
  
   end 
% 
%     if band ==1 % CA
%          TOARad_L8_all.Band1=MeanTOArad_L8';
%      elseif band==2 % Blue
%          TOARad_L8_all.Band2=MeanTOArad_L8';
%      elseif band==3 % Green
%          TOARad_L8_all.Band3=MeanTOArad_L8';
%      elseif band==4 % Red
%          TOARad_L8_all.Band4=MeanTOArad_L8';
%      elseif band==5 % NIR
%          TOARad_L8_all.Band5=MeanTOArad_L8';
%      elseif band==6 % SWIR1
%          TOARad_L8_all.Band6=MeanTOArad_L8';
%      elseif band==7 % SWIR2
%          TOARad_L8_all.Band7=MeanTOArad_L8';
%    end
   
   end 