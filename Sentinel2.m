%%
clc, close, clear all;

base= 'Z:\ImageDrive\Sentinel\MSI-A\P181\R040';
dates = dir(base);
dates([1 2])=[];

%ROI coordinates
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

Q_value=10000;
bands=  {'01' '02' '03' '04' '8A' '11' '12'}; 
band_name= {'CA','B',  'G', 'R', 'NIR', 'SWIR1', 'SWIR2'};

for band= 2%:7
 for date = 2 %1:size(dates,1)
   if exist (fullfile(base, dates(date).name, 'MSIL1C'),'dir')
        
        %path of the xml file
        MSIpath = dir(fullfile(base,dates(date).name,'MSIL1C','*MSIL1C.xml'));
        MTLpath = dir(fullfile(base,dates(date).name,'MSIL1C','*TL.xml'));

        MSIfile=fullfile(base, dates(date).name,'MSIL1C', MSIpath.name);
        MTDfile=fullfile(base, dates(date).name,'MSIL1C', MTLpath.name);

        %Storing the values from xml file
        [MSI_values]= xml2struct_new_v(MSIfile);
        [MTD_TL] = xml2struct_new_v(MTDfile);
        
        % Image directory and Image file name
        Band_No(band)=dir(fullfile(base,dates(date).name,'MSIL1C', strcat('*B', bands{band},'.jp2')));
        S2_image_file= strcat(base, filesep, dates(date).name, filesep, 'MSIL1C', filesep, Band_No(band).name);
        
        %BandID selection depending on band number
        if band == 2 || band ==3 ||band == 4 || band == 8
            BandID = 1; %% 10m

        elseif band == 1 ||band == 9 ||band == 10
            BandID = 3; %% 60m

        elseif  band == 5 || band == 6 ||band == 7
            BandID = 2; %% 20m
        end

        % Extrating Image infos depending on band number which is defined by BandID
        s2_image_info.CornerCoords.X=str2num(MTD_TL.n1_colon_Level_dash_1C_Tile_ID....
                    .n1_colon_Geometric_Info.Tile_Geocoding.Geoposition{1, BandID}.ULX.Text);
                
        s2_image_info.CornerCoords.Y=str2num(MTD_TL.n1_colon_Level_dash_1C_Tile_ID....
                    .n1_colon_Geometric_Info.Tile_Geocoding.Geoposition{1, BandID}.ULY .Text);

        s2_image_info.Height=str2num(MTD_TL.n1_colon_Level_dash_1C_Tile_ID....
                    .n1_colon_Geometric_Info.Tile_Geocoding.Size{1, BandID}.NROWS.Text);

        s2_image_info.Width=str2num(MTD_TL.n1_colon_Level_dash_1C_Tile_ID....
                    .n1_colon_Geometric_Info.Tile_Geocoding.Size{1, BandID}.NCOLS.Text);

        s2_image_info.CurrentRes=str2num(MTD_TL.n1_colon_Level_dash_1C_Tile_ID.n1_colon_Geometric_Info....
                    .Tile_Geocoding.Geoposition{1, BandID}.Attributes.resolution);
                    
        %Reference Matrix 
        R_S2.RefMatrix=[0 -s2_image_info.CurrentRes; s2_image_info.CurrentRes 0; s2_image_info.CornerCoords.X  s2_image_info.CornerCoords.Y];
        
        %Map coordinate to pixel coordinate
        [Pixel_Row_unrounded, Pixel_Column_unrounded] = map2pix(R.RefMatrix, x_vec, y_vec);
        Pixel_Row= round(Pixel_Row_unrounded);
        Pixel_Column= round(Pixel_Column_unrounded);
        
        mask= poly2mask(Pixel_Column, Pixel_Row, s2_image_info.Height, s2_image_info.Width);
        
        %Image Reading
        Image = imread(S2_image_file);   
        Image = double(Image);
        
        %TOA Reflectance
        Image_mask=Image.*mask;
        TOARef=(Image_mask)/Q_value;
        TOARef(TOARef==0)= nan;  
        
%         TOARef = TOARef(~isnan(TOARef));
%         meanTOARef(date)=mean2(TOARef(~isnan(TOARef)));
        
        % solar irradiance is given as bandId 0 to 12 ; BandId 0 represents band1 (from MSI file)
        % here bandId="0" is CA band, bandId="7" is Band 8-NIR band,
        % bandID="8" is Band 8A- Narrow NIR % bandID="9" is Band 9- -water
        % vapor same onwards
       
        if band==1
            bandID=1; %CA band
        elseif band==2
            bandID=2; %Blue band
        elseif band==3
            bandID=3; %Green band
        elseif band==4
            bandID=4; %Red Band
        elseif band==5 % for band8A (comparable to L8-NIR)
            bandID=9; % 
        elseif band==6 %for band 11 (comparable to L8-SWIR1)
            bandID=12; % don't get confused the xml files bandID is just an index
        elseif band==7 % for band 12 (comparable to L8-SWIR2)
            bandID=13;
        end
        
        %Solar Irradiance
        TOArad_con_parameter.ESUN=str2double(MSI_values.n1_colon_Level_dash_1C_User_Product.n1_colon_General_Info.....
             .Product_Image_Characteristics.Reflectance_Conversion.Solar_Irradiance_List.SOLAR_IRRADIANCE{bandID}.Text);
        TOArad_con_parameter.d=str2double(MSI_values.n1_colon_Level_dash_1C_User_Product.n1_colon_General_Info.....
             .Product_Image_Characteristics.Reflectance_Conversion.U.Text);
         
        %without considering earth-sun distance
        %TOARadS2b1(date)=TOARef*TOArad_con_parameter.ESUN/(pi);
        
        %TOArad_S2=TOARef*TOArad_con_parameter.ESUN/(TOArad_con_parameter.d.^2*pi);
       
        TOArad_S2=TOARef*TOArad_con_parameter.ESUN/(pi);
        TOArad_S2_sampled=imresize(TOArad_S2, 10/30); % downsampling
         
        %MeanTOArad_S2(date)=meanTOARef(date)*TOArad_con_parameter.ESUN/(pi);
        % Day of year and Decimal year
        acquisition_date =MTD_TL.n1_colon_Level_dash_1C_Tile_ID....
                    .n1_colon_General_Info.SENSING_TIME.Text(end-23:end-14);
        [doy,fraction] = date2doy(datenum(acquisition_date));
        DoY(date)=doy; %day of year
        DateVec= datevec(acquisition_date);
        DeciYear(date)=DateVec(1,1)+DateVec(1,2)./12+DateVec(1,3)./365; %decimal year   
      
  end
       
 end
 
%     if band ==1 % CA
%          TOARad_S2_all.Band1=MeanTOArad_S2';
%      elseif band==2 % Blue
%          TOARad_S2_all.Band2=MeanTOArad_S2';
%      elseif band==3 % Green
%          TOARad_S2_all.Band3=MeanTOArad_S2';
%      elseif band==4 % Red
%          TOARad_S2_all.Band4=MeanTOArad_S2';
%      elseif band==5 % NIR
%          TOARad_S2_all.Band5=MeanTOArad_S2';
%      elseif band==6 % SWIR1
%          TOARad_S2_all.Band6=MeanTOArad_S2';
%      elseif band==7 % SWIR2
%          TOARad_S2_all.Band7=MeanTOArad_S2';
%     end
   
end

%%
%not using

new_R.RefMatrix(1,2)=-30;
new_R.RefMatrix(2,1)=30;

load('l8_info.mat')

load('L8_R.mat')
load('TOArad_L8.mat')
load('TOAref_L8.mat') % for checking only

%%%
[nR,nC]=map2pix(R, s2_image_info.CornerCoords.X, s2_image_info.CornerCoords.Y);

TOArad_S2_B2 = NaN(size(TOArad_L8,1), size(TOArad_L8,2));
TOArad_S2_B2(nR:nR+3660-1,nC:nC+3660-1)= TOArad_S2_sampled;



TOArad_L8(nR:nR+3660-1,nC:nC+3660-1)=t;
TOArad_S2_c=TOArad_L8()
TOArad_S2_c(nR:nR+3660-1,nC:nC+3660-1)=t;

%%%checking
S2_TOAref_b6 = NaN(size(TOAref_L8,1), size(TOAref_L8,2));
S2_TOAref_b6(nR:nR+3660-1,nC:nC+3660-1)=tref;


% Load_S2_check=load('S2_Data.mat');
% TOARad_S2_Final.all= [Load_S2_check.TOARad_S2_all.Band1 Load_S2_check.TOARad_S2_all.Band2 Load_S2_check.TOARad_S2_all.Band3....
%                     Load_S2_check.TOARad_S2_all.Band4 Load_S2_check.TOARad_S2_all.Band5 Load_S2_check.TOARad_S2_all.Band6....
%                     Load_S2_check.TOARad_S2_all.Band7];
% 
% Deci= find(isnan(TOARad_S2_Final.all(:,1)));
% DeciYear_Final.all= Load_S2_check.DeciYear'; 
% 
% TOARad_S2_all.Band1  
% TOARad_S2_Final.all(~any(~isnan(TOARad_S2_Final.all), 2), :)=[];
% 
% % TOARad_S2_Final2.NaN(6,:)
% % 
% % DeciYear_Final.alldates=DeciYear';
% 
% for D = 1 : length(Deci)
%    De= Deci(D);
%    DeciYear_Final.all(De)=0;
% end
% 
% DeciYear_Final.all(find(~DeciYear_Final.all))=[];
% 
% for b = 1:7
%     outlier = abs(TOARad_S2_Final.all(:,b) - mean(TOARad_S2_Final.all(:,b))) > 2*std(TOARad_S2_Final.all(:,b));
%     outlier_index = find(outlier);
% end
% 
% outlier1 = abs(TOARad_S2_Final.all(:,1))-mean(TOARad_S2_Final.all(:,1))>2*std(TOARad_S2_Final.all(:,1));
% outlier_index1 = find(outlier1);
% TOARad_S2_Final.all(outlier_index,:)=[]

% 
% for band= 2:7
%     TOARad_S2.NaN(isnan(TOARad_S2.NaN))=[];  
% end

% % storing all the TOA Radiances in a table
% Table=struct2table(TOARad_S2_all); % taking all the values to table (type table)
% vars = {'Band1' 'Band2' 'Band3' 'Band4' 'Band5' 'Band6' 'Band7'};
% Table2= Table{:,vars}; % changing table type to double for isnan
% %N2cell = num2cell(Table2); % can be done either way
% T2Cell = table2cell(Table); %same as previous line
% T2Cell(isnan(Table2)) ={'NaN'}; % to export all the NaN values to excel % isnan doesn't take table type
% xlswrite('TOARad_S2.xlsx', T2Cell);
% 
% mytable = rmmissing(A2); 
