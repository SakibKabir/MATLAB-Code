%UPDATE: 22 AUGUST, NAHID; WILL TRY TO STAMP THE NEW TILE BASED COOKIE
%CUTTERS %woops run ion windows please. unix is giving error while reading
%the xml file

%UPDATE: 30 JULY,2018: DEVELOP THE LOGIC OF 3 DIFFERENT RESOLUTION COOKIE
%CUTTER MASK

%UPDATE:27 JULY,2018. PRESERVING PATH ROW STRUCTURE SIMILLAR TO LANDSAT
%SERIES 

%PREVIOUS UPDATES: DATE: LOST
% WORK ON PROGRESS
% SPECIFY THE PATH AND ROW AS PAIR BY PAIR TO BE PROCESSED
% CLUSTER13 SELECTED PATH ROWS ARE BELLOW FOR LANDSAT 8
% CLSUTER13:[178 47;179 41;180 40;181 40;182 40;185 47;186 47;188 47;
% 187 47;189 46;190 43;191 37;192 37;193 37;199 46;200 47]

clc
clear all
close all
tic

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

%PREVIOUS PATH ROW LIST 
path_row_pairs=[181 40];


%SCALE FACTOR FOR TOA REFLECTANCE
Q_value=10000;
angle_factor=10;
resolution=[60   10    10   10    20    20      20];
bands=    {'01' '02'  '03' '04' '8A'  '11'     '12'}; %literally would be b
band_str= {'CA','B',  'G', 'R', 'NIR', 'SWIR1', 'SWIR2'};
%IMP: MENTION THE COOKIE CUTTER FILE BELLOW IF NOT CHOSEN AUTOMATICALLY 
%%%%%%%%%%%%%
num_path_rows=size(path_row_pairs,1);

date=1
base= 'Z:\ImageDrive\Sentinel\MSI-A\P181\R040';
dates = dir(base);
dates([1 2])=[];
        MSIpath = dir(fullfile(base,dates(date).name,'MSIL1C','*MSIL1C.xml'));
        MTDpath = dir(fullfile(base,dates(date).name,'MSIL1C','*TL.xml'));
        
        MSIfile=fullfile(base, dates(date).name,'MSIL1C', MSIpath.name);
        MTDfile=fullfile(base, dates(date).name,'MSIL1C', MTDpath.name);
        [xml_values] = xml2struct_new_v(MSIfile);
        [MTD_TL] = xml2struct_new_v(MTDfile);

for path_row_count=1:num_path_rows%start of path_row_pair loop

    clearvars -except path_row_pairs path_row_count num_path_rows bands resolution band_str Q_value angle_factor;   
    %%%%%%%%%create path_row directory of target landsat images
    if ispc
        image_drive_location = strcat('Z:\ImageDrive\Sentinel\MSI-A\P'...
            ,sprintf('%03d',path_row_pairs(path_row_count,1)),'\R'...
            ,sprintf('%03d',path_row_pairs(path_row_count,2)));
    elseif isunix  
        image_drive_location = strcat('~/zdrive/ImageDrive/Sentinel/MSI-A/P'...
            ,sprintf('%03d',path_row_pairs(path_row_count,1)),'/R'...
            ,sprintf('%03d',path_row_pairs(path_row_count,2)));   
    end
    %%%%%%%%%%
    %make the cookie cutter address strings for 3 differnt resoulutions
    temp_counter=0;
    for res=[10 20 60]
         temp_counter=temp_counter+1;
         cookie_cutter_file{temp_counter}=strcat('cookie_cutters',filesep,'SaharanImageDisplay19Cluster13_CC_Local_'...
        ,sprintf('%03d',path_row_pairs(path_row_count,1))...
        ,sprintf('%02d',path_row_pairs(path_row_count,2)),'_',num2str(res),'.tif');
    end
    clear temp_counter res

    acquisition_dates_available = dir(image_drive_location);

    %DELETE EMPTY STRUCTURES
    acquisition_dates_available([1,2])=[];
    num_of_acquisition_dates= length(acquisition_dates_available(:,1)); 

    %%%%folder to save staistics with angle information
    % % save_path=strcat('ruchiras_stat_path',sprintf('%03d',path_row_pairs(path_row_count,1)),'_','row',sprintf('%03d',path_row_pairs(path_row_count,2)));
    % % mkdir(save_path);
    %%%%%%

%%%%%%%%start processing the dates in a path and row folder
    for i = 1%:num_of_acquisition_dates
        
        date= strcat(image_drive_location,filesep,acquisition_dates_available(i).name);%creating the date directory string
        if exist(date, 'file')
        %disp('!! Date Folder exists')
        else
            continue
        end
        ImageLocation = strcat(date,filesep,'MSIL1C',filesep); %create path to LC1 folder
        %check whether the LC1 folder exists or not
        if ( exist(ImageLocation, 'dir') ~= 7)
            display('L1C folder doesnt exists !!! Moving to next acquisition !!!');
            continue

        end

        ListOfFilesOnThisDate = dir(ImageLocation);
        %get all lc1 folders from all dates available
        
          if size(ListOfFilesOnThisDate,1)<3 %%%checking if the LC1 folder is empty
            continue
          end
          ListOfFilesOnThisDate([1 2])=[];
        
          
        %%%%finding base file name
        %%%%will target the jp2 extension
        for j=1:size(ListOfFilesOnThisDate,1)
            base_flag=strcmp(ListOfFilesOnThisDate(j).name(end-2:end),'jp2');
            if base_flag==1
                break
            end
        end
        
        BasenameS2AB = ListOfFilesOnThisDate(j).name(1:end-8);%extracting comoon basename for the data 
        
        MTLFilenameL8 = fullfile(ImageLocation(1,:),'MTD_TL.xml');
        
          if ~(exist(MTLFilenameL8))
            continue
          end
          
        MTLList = xml2struct_new_v(MTLFilenameL8);%read the mtl file

        
        %read sun and sensor angles 
        try
            SolarZenith = imread(strcat(ImageLocation(1,:),'Sun_Zenith.tif'));
            SolarAzimuth= imread(strcat(ImageLocation(1,:),'Sun_Azimuth.tif'));

            SensorZenith =imread(strcat(ImageLocation(1,:),'View_Zenith.tif'));
            SensorAzimuth =imread(strcat(ImageLocation(1,:),'View_Azimuth.tif'));
        catch
            continue
        end

%%give a progress indicator 
disp(['Paths/Rows (' num2str((path_row_count/num_path_rows)*100) '%): '    num2str(path_row_count) '/' num2str(num_path_rows)])
disp(['  Dates being processed (' num2str((i/num_of_acquisition_dates)*100) '%): ' num2str(i) '/' num2str(num_of_acquisition_dates)])


           %BEGIN THREE RESOLUTION PROCESSING
           for c=[1 2 3] %1 FOR 10 METER; 2 FOR 20 METER  AND 3 FOR 60 METER ;
                s2_image_info.CornerCoords.X=str2num(MTLList.n1_colon_Level_dash_1C_Tile_ID....
                    .n1_colon_Geometric_Info.Tile_Geocoding.Geoposition{1, c}.ULX.Text);

                s2_image_info.CornerCoords.Y=str2num(MTLList.n1_colon_Level_dash_1C_Tile_ID....
                    .n1_colon_Geometric_Info.Tile_Geocoding.Geoposition{1, c}.ULY .Text);

                s2_image_info.Height=str2num(MTLList.n1_colon_Level_dash_1C_Tile_ID....
                    .n1_colon_Geometric_Info.Tile_Geocoding.Size{1, c}.NROWS.Text);

                s2_image_info.Width=str2num(MTLList.n1_colon_Level_dash_1C_Tile_ID....
                    .n1_colon_Geometric_Info.Tile_Geocoding.Size{1, c}.NCOLS.Text);

                s2_image_info.CurrentRes=str2num(MTLList.n1_colon_Level_dash_1C_Tile_ID.n1_colon_Geometric_Info....
                    .Tile_Geocoding.Geoposition{1, c}.Attributes.resolution);

           
            switch c

               case 1
                    cookie_cutter=geotiffread(cookie_cutter_file{c});
                    cookie_cutter_info=geotiffinfo(cookie_cutter_file{c});

                   [Pixel_Row_unrounded, Pixel_Col_unrounded] = map2pix(cookie_cutter_info.RefMatrix,s2_image_info.CornerCoords.X,s2_image_info.CornerCoords.Y);%ls_image_info.RefMatrix(3,1),ls_image_info.RefMatrix(3,2));
                   Pixel_Row=round(Pixel_Row_unrounded); %pixel column/row number can not be fraction
                   Pixel_Col=round(Pixel_Col_unrounded);
                   try
                   cookie_cutter_roi=cookie_cutter(Pixel_Row(1,1):Pixel_Row(1,1)+(s2_image_info.Height-1),Pixel_Col(1,1):Pixel_Col(1,1)+(s2_image_info.Width-1));%crop of original cookie_cutter
                   
                   catch
                       continue
                   end
                   cookie_cutter_roi=double(cookie_cutter_roi);
                  %%display progress
                   clc
                     disp(['Paths/Rows (' num2str((path_row_count/num_path_rows)*100) '%): '    num2str(path_row_count) '/' num2str(num_path_rows)])
                     disp(['  Dates being processed (' num2str((i/num_of_acquisition_dates)*100) '%): ' num2str(i) '/' num2str(num_of_acquisition_dates)])
                  % disp(['  Bands being processed: ' num2str(k) '/' '7'])%progress indicator
                  %%display progress
                    for k= [2 3 4] %THESE ARE 10M BAND
                        S2_image_name = strcat(ImageLocation(1,:),BasenameS2AB,'_B', bands{k},'.jp2');%band by band name creation
                        ImageIn = imread(S2_image_name);     
                        ImageIn = double(ImageIn);          
                        ImageIn(~ImageIn) = nan;%Convert 0s to NaN
                        TOARef=(ImageIn.*cookie_cutter_roi)/Q_value;
                        TOARef(TOARef==0) = nan; 
                        meanTOA(i,k) = mean2(TOARef(~isnan(TOARef)));
                        stdTOA(i,k) = std2(TOARef(~isnan(TOARef)));
                        
                        
                        if k==2 %LIMIT THE ANGLE  PROCESSING TO ONLY ONE BAND 
                            SolarZenith=(double(SolarZenith).*cookie_cutter_roi)/angle_factor;
                            SolarAzimuth=(double(SolarAzimuth).*cookie_cutter_roi)/angle_factor;
                            SensorZenith=(double(SensorZenith).*cookie_cutter_roi)/angle_factor;
                            SensorAzimuth=(double(SensorAzimuth).*cookie_cutter_roi)/angle_factor;
                        end
                        
                    end

                
                
                case 2
                    
                    cookie_cutter=geotiffread(cookie_cutter_file{c});
                    cookie_cutter_info=geotiffinfo(cookie_cutter_file{c});
                    [Pixel_Row_unrounded, Pixel_Col_unrounded] = map2pix(cookie_cutter_info.RefMatrix,s2_image_info.CornerCoords.X,s2_image_info.CornerCoords.Y);%ls_image_info.RefMatrix(3,1),ls_image_info.RefMatrix(3,2));
                    Pixel_Row=round(Pixel_Row_unrounded); %pixel column/row number can not be fraction
                    Pixel_Col=round(Pixel_Col_unrounded);
                    
                    try
                    cookie_cutter_roi=cookie_cutter(Pixel_Row(1,1):Pixel_Row(1,1)+(s2_image_info.Height-1),Pixel_Col(1,1):Pixel_Col(1,1)+(s2_image_info.Width-1));%crop of original cookie_cutter
                    catch
                       continue
                    end
                    
                    cookie_cutter_roi=double(cookie_cutter_roi);    
                    
                    
                    
                    clc
                    for k= [5 6 7] %THESE ARE 20M BAND 
                        S2_image_name = strcat(ImageLocation(1,:),BasenameS2AB,'_B', bands{k},'.jp2');%band by band name creation
                        ImageIn = imread(S2_image_name);     
                        ImageIn = double(ImageIn);          
                        ImageIn(~ImageIn) = nan;%Convert 0s to NaN
                        TOARef=(ImageIn.*cookie_cutter_roi)/Q_value;
                        TOARef(TOARef==0) = nan; 
                        meanTOA(i,k) = mean2(TOARef(~isnan(TOARef)));
                        stdTOA(i,k) = std2(TOARef(~isnan(TOARef)));
                    end
                
                case 3 %STANDS FOR 6O METER
                    
                   %cookie_cutter=geotiffread(cookie_cutter_file{c});
                   %cookie_cutter_info=geotiffinfo(cookie_cutter_file{c});
                   cookie_cutter_info.RefMatrix=[0 -60;60 0;s2_image_info.CornerCoords.X  s2_image_info.CornerCoords.Y];
                   [Pixel_Row_unrounded, Pixel_Col_unrounded] = map2pix(cookie_cutter_info.RefMatrix,x_vec,y_vec);
                   Pixel_Row=round(Pixel_Row_unrounded); %pixel column/row number can not be fraction
                   Pixel_Col=round(Pixel_Col_unrounded);
                   mask60=poly2mask(Pixel_Row,Pixel_Col,s2_image_info.Height,s2_image_info.Width); 
                   try
                   %cookie_cutter_roi=cookie_cutter(Pixel_Row(1,1):Pixel_Row(1,1)+(s2_image_info.Height-1),Pixel_Col(1,1):Pixel_Col(1,1)+(s2_image_info.Width-1));%crop of original cookie_cutter
                   catch
                       continue
                   end
                   cookie_cutter_roi=double(cookie_cutter_roi);  
                   
                    clc
                    for k= [1] %IT IS 60M BAND
                        S2_image_name = strcat(ImageLocation(1,:),BasenameS2AB,'_B', bands{k},'.jp2');%band by band name creation
                        
                        ImageIn = imread(S2_image_name);     
                        ImageIn = double(ImageIn);          
                        ImageIn(~ImageIn) = nan;%Convert 0s to NaN
                        TOARef=(ImageIn.*mask60)/Q_value;
                        TOARef(TOARef==0) = nan;
                       
                        meanTOAref = mean2(TOARef(~isnan(TOARef)));
                        stdTOA(i,k) = std2(TOARef(~isnan(TOARef)));
                         %%%%%%
                        bandID=k;
                         x=xml_values.n1_colon_Level_dash_1C_User_Product.n1_colon_General_Info.Product_Image_Characteristics.Reflectance_Conversion.Solar_Irradiance_List.SOLAR_IRRADIANCE;
                         ESUN=str2double(x{1,bandID}.Text);

                         dt=str2double(xml_values.n1_colon_Level_dash_1C_User_Product.n1_colon_General_Info.Product_Image_Characteristics.Reflectance_Conversion.U.Text);
                         d=sqrt(1/dt);
                         MeanTOArad=meanTOAref*ESUN/(pi*d.^2);
                        %%%%%%
                    end
           end
          %SWITCH ENDS HERE 
          
          
           end
          %END THREE RESOLUTION PROCESSING HERE
           

                error_counter=0;
                try
                stat(i).mean=meanTOA(i,:);
                stat(i).std=stdTOA(i,:);
                stat(i).spat_unc=100*(stdTOA(i,:)./meanTOA(i,:));
                stat(i).acq_date=acquisition_dates_available(i).name;
                catch
                error_counter=error_counter+1;
                errorflag{error_counter}=strcat('tiles might be mismatched at path',sprintf('%03d',path_row_pairs(path_row_count,1)),'_','row',sprintf('%03d',path_row_pairs(path_row_count,2)));
                continue    
                end
                
                
                %%%%% embedd the mean values in data file
                SolarZenith(SolarZenith==0)=nan;
                SolarAzimuth(SolarAzimuth==0)=nan;
                SensorZenith(SensorZenith==0)=nan;
                SensorAzimuth(SensorAzimuth==0)=nan;
                
                %%%%%
                
                %%%embedd them in stat structure too
                stat(i).SolZenMean= mean2(SolarZenith(~isnan(SolarZenith)));
                stat(i).SolAzmMean=mean2(SolarAzimuth(~isnan(SolarAzimuth)));
                stat(i).SenZenMean= mean2(SensorZenith(~isnan(SensorZenith)));
                stat(i).SenAzmMean= mean2(SensorAzimuth(~isnan(SensorAzimuth)));
                

% % %   save(strcat(save_path,filesep,'cluster13_stat_path',sprintf('%03d',path_row_pairs(path_row_count,1)),...
% % %   '_','row',sprintf('%03d',path_row_pairs(path_row_count,2)),'_',...
% % %   sprintf('%04d',i)),'stat_with_angle','-v7.3');
% % %   clear stat_with_angle
               
                
end 
%%%%%%%end processing the dates in a path and row folder
%stat(1:2)=[];
save(strcat('August24_sentinel_2B_c13_stat_path',sprintf('%03d',path_row_pairs(path_row_count,1)),'_','row',sprintf('%03d',path_row_pairs(path_row_count,2))),'stat','-v7.3');%later include path row info

end %end of path_row pair loop
disp('done')
toc
%exit

