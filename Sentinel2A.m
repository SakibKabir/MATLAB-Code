clc, close, clear all;

base= 'Z:\ImageDrive\Sentinel\MSI-A\P181\R040';
dates = dir(base);

UL_y=3162691;
UL_x=739080;
UR_y=3162576;
UR_x=770212;
LR_y=3126914;
LR_x=770498;
LL_y=3128004;
LL_x=731970;

x_vec=[UL_x UR_x LR_x LL_x UL_x];
y_vec=[UL_y UR_y LR_y LL_y UL_y];

Resolution=[60 10 10 10 20 20 20 10 60 60 20 20 20];
Q_value=10000;
Refl_S2A=zeros(length(1),7);

      
dates([1 2])=[];
%    for date = 1:size(dates,1)
   for date = 1
   if exist (fullfile(base, dates(date).name, 'MSIL1C'),'dir')
       
       
        Band1=dir(fullfile(base,dates(date).name,'MSIL1C','*B01.jp2'));
        
        S2_image_name=fullfile(base,dates(date).name,'MSIL1C',Band1.name);
        
        %path of the xml file
        MSIpath = dir(fullfile(base,dates(date).name,'MSIL1C','*MSIL1C.xml'));
        MTLpath = dir(fullfile(base,dates(date).name,'MSIL1C','*TL.xml'));
        
        MSIfile=fullfile(base, dates(date).name,'MSIL1C', MSIpath.name);
        MTDfile=fullfile(base, dates(date).name,'MSIL1C', MTLpath.name);
        
        %Storing the values from xml file
        [MSI_values]= xml2struct_new_v(MSIfile);
        [MTD_TL] = xml2struct_new_v(MTDfile);
        
        
%         % Day of Year and Decimal Year
%         acquisition_date =MTL_List.PRODUCT_METADATA.DATE_ACQUIRED;
%         [doy,fraction] = date2doy(datenum(acquisition_date));
%         DoY(date)=doy; %day of year
%         DateVec= datevec(acquisition_date);
%         DeciYear(date)=DateVec(1,1)+DateVec(1,2)./12+DateVec(1,3)./365; %decimal year

        %Extracting the necsseary values from xml file
         bandID=3; 
         s2_image_info.CornerCoords.X=str2num(MTD_TL.n1_colon_Level_dash_1C_Tile_ID....
                    .n1_colon_Geometric_Info.Tile_Geocoding.Geoposition{1, bandID}.ULX.Text);
                
         s2_image_info.CornerCoords.Y=str2num(MTD_TL.n1_colon_Level_dash_1C_Tile_ID....
                    .n1_colon_Geometric_Info.Tile_Geocoding.Geoposition{1, bandID}.ULY .Text);

         s2_image_info.Height=str2num(MTD_TL.n1_colon_Level_dash_1C_Tile_ID....
                    .n1_colon_Geometric_Info.Tile_Geocoding.Size{1, bandID}.NROWS.Text);

         s2_image_info.Width=str2num(MTD_TL.n1_colon_Level_dash_1C_Tile_ID....
                    .n1_colon_Geometric_Info.Tile_Geocoding.Size{1, bandID}.NCOLS.Text);

         s2_image_info.CurrentRes=str2num(MTD_TL.n1_colon_Level_dash_1C_Tile_ID.n1_colon_Geometric_Info....
                    .Tile_Geocoding.Geoposition{1, bandID}.Attributes.resolution);
                
        %Reference Matrix 
        R.RefMatrix=[0 -s2_image_info.CurrentRes; s2_image_info.CurrentRes 0; s2_image_info.CornerCoords.X  s2_image_info.CornerCoords.Y];
        
        %Map coordinate to pixel coordinate
        [Pixel_Row_unrounded, Pixel_Column_unrounded] = map2pix(R.RefMatrix, x_vec, y_vec);
        Pixel_Row= round(Pixel_Row_unrounded);
        Pixel_Column= round(Pixel_Column_unrounded);
        
        mask= poly2mask(Pixel_Row, Pixel_Column, s2_image_info.Height, s2_image_info.Width);
        
        Image = imread(S2_image_name);     
        Image = double(Image);          
        Image(~Image) = nan;  %Convert 0s to NaN
        TOARef=(Image.*mask)/Q_value;
        TOARef(TOARef==0)= nan;  
        TOARef = TOARef(~isnan(TOARef));
        meanTOARef=mean2(TOARef(~isnan(TOARef)));
         
         band=1; % solar irradiance is given as BandId 0 to 12 ; BandId 0 represents band1
         TOArad_con_parameter.ESUN=str2double(MSI_values.n1_colon_Level_dash_1C_User_Product.n1_colon_General_Info.....
             .Product_Image_Characteristics.Reflectance_Conversion.Solar_Irradiance_List.SOLAR_IRRADIANCE{1,bandID}.Text);
         TOArad_con_parameter.d=str2double(MSI_values.n1_colon_Level_dash_1C_User_Product.n1_colon_General_Info.....
             .Product_Image_Characteristics.Reflectance_Conversion.U.Text);
         
         %without considering earth-sun distance
         TOARadS2b1=TOARef*TOArad_con_parameter.ESUN/(pi);
         MeanTOArad=meanTOARef*TOArad_con_parameter.ESUN/(pi);
         
         %d=sqrt(1/TOArad_con_parameter.dt);
         %MeanTOArad=meanTOARef*TOArad_con_parameter.ESUN/(pi*TOArad_con_parameter.d.^2);
         
         
%        ToA4Mean(date) = mean(ToA4ROI(ToA4ROI~=0));
%        Sun_Azimuth(date)= MTL_List.GROUP_IMAGE_ATTRIBUTES.SUN_AZIMUTH;
%        Sun_Elevation(date)=MTL_List.GROUP_IMAGE_ATTRIBUTES.SUN_ELEVATION;
%        ToA4Mean_corr(date)=ToA4Mean(date)./sind(Sun_Elevation(date));
       
 
   end

   end


 