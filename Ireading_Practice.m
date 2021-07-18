close, clear all;

base= 'Z:\ImageDrive\OLI.TIRS\L8\P181\R040';
dates = dir(base);
dates([1 2])=[];

  for date = 1:size(dates,1)
      if exist (fullfile(base, dates(date).name, 'LC1'),'dir')
        Band1=dir(fullfile(base,dates(date).name,'LC1','*B1.tif')); % deep blue and violets
        Band2=dir(fullfile(base,dates(date).name,'LC1','*B2.tif')); % visible blue
        Band3=dir(fullfile(base,dates(date).name,'LC1','*B3.tif')); % green
        Band4=dir(fullfile(base,dates(date).name,'LC1','*B4.tif')); % red
        Band5=dir(fullfile(base,dates(date).name,'LC1','*B5.tif')); % NIR
        Band6=dir(fullfile(base,dates(date).name,'LC1','*B6.tif')); % SWIR1
        Band7=dir(fullfile(base,dates(date).name,'LC1','*B7.tif')); % SWIR2
        Band8=dir(fullfile(base,dates(date).name,'LC1','*B8.tif')); % Pan
        Band9=dir(fullfile(base,dates(date).name,'LC1','*B9.tif')); % Cirrus
        
%       Band10=dir(fullfile(base,dates(date).name,'LC1','*B10.tif'));% TIR
%       Band11=dir(fullfile(base,dates(date).name,'LC1','*B11.tif'));% TIR
         
        % MTL Parser Function to extract all the data from the Meta Data File
        MTL=dir(fullfile(base, dates(date).name,'LC1','*MTL.txt'));
        [MTL_List, vaule]= MTL_parser_L8(fullfile(MTL.folder, MTL.name));
        
        % Day of Year and Decimal Year
        acquisition_date =MTL_List.PRODUCT_METADATA.DATE_ACQUIRED;      %Date extraction from MetaData file
        [doy,fraction] = date2doy(datenum(acquisition_date));           %date2doy function to convert date to get Day of the Year
        DoY(date)=doy;                                                  %day of year
        DateVec= datevec(acquisition_date);                             %converting acquisition date to date vector
        DeciYear(date)=DateVec(1,1)+DateVec(1,2)./12+DateVec(1,3)./365; %decimal year
        
        % Interested Regions Lat and Lon Information
        UL_LongitudeL8=727396.6113;
        UL_LatitudeL8=3234770.0903;
        UR_LongitudeL8=810868.7842;
        UR_LatitudeL8=3234207.0337;
        LR_LongitudeL8=812276.4111;
        LR_LatitudeL8=3168330.1758;
        LL_LongitudeL8=726974.3188;
        LL_LatitudeL8=3168189.4189;
        
        %----------------Band1:costal/aerosol band------------------------%
        rmb1(date)= MTL_List.RADIOMETRIC_RESCALING.RADIANCE_MULT_BAND_1;
        rab1(date)= MTL_List.RADIOMETRIC_RESCALING.RADIANCE_ADD_BAND_1;
        [A1, R1] = geotiffread(fullfile(base,dates(date).name,'LC1',Band1.name));
        A1=double(A1);
        %DN1(date)=mean2(A1(r,c));
        ToA_rad1=A1*rmb1(date)+rab1(date);  % ToA Radiance calculation                       
      
        % Calculation of pixel coordinates from map coordinates
        [UL_Row, UL_Column] = map2pix(R1, UL_LongitudeL8, UL_LatitudeL8);
        [UR_Row, UR_Column] = map2pix(R1, UR_LongitudeL8, UR_LatitudeL8);
        [LL_Row, LL_Column] = map2pix(R1, LL_LongitudeL8, LL_LatitudeL8);
        [LR_Row, LR_Column] = map2pix(R1, LR_LongitudeL8, LR_LatitudeL8);
        
        ROIRow1 = [UL_Row UR_Row LR_Row LL_Row UL_Row]; 
        ROIColumn1 = [UL_Column UR_Column LR_Column LL_Column UL_Column];
        [Row, Column]= size(A1);
        
        % Region of Interest(ROI) Mask 
        mask= poly2mask(ROIColumn1, ROIRow1, Row, Column);
        ToA1ROI = double(ToA_rad1).*mask;
        
        %Mean of the ToA Radiance of Interested Region
        ToA1Mean(date) = mean(ToA1ROI(ToA1ROI~=0));
        
        %ToA Radiance with a correction for Sun Angle
        Sun_Azimuth(date)= MTL_List.GROUP_IMAGE_ATTRIBUTES.SUN_AZIMUTH;
        Sun_Elevation(date)=MTL_List.GROUP_IMAGE_ATTRIBUTES.SUN_ELEVATION;
        ToA1Mean_corr(date)=ToA1Mean(date)./sind(Sun_Elevation(date));
        
        %-----------------------Band2:Visible blue------------------------%
        rmb2(date)= MTL_List.RADIOMETRIC_RESCALING.RADIANCE_MULT_BAND_2;
        rab2(date)= MTL_List.RADIOMETRIC_RESCALING.RADIANCE_ADD_BAND_2;
        [A2, R2] = geotiffread(fullfile(base,dates(date).name,'LC1',Band2.name));
        A2=double(A2);
        %DN2(date)=mean2(A2(r,c));
        ToA_rad2=A2*rmb2(date)+rab2(date);  % ToA Radiance calculation                       
      
        % Calculation of pixel coordinates from map coordinates
        [UL_Row, UL_Column] = map2pix(R2, UL_LongitudeL8, UL_LatitudeL8);
        [UR_Row, UR_Column] = map2pix(R2, UR_LongitudeL8, UR_LatitudeL8);
        [LL_Row, LL_Column] = map2pix(R2, LL_LongitudeL8, LL_LatitudeL8);
        [LR_Row, LR_Column] = map2pix(R2, LR_LongitudeL8, LR_LatitudeL8);
        
        ROIRow2 = [UL_Row UR_Row LR_Row LL_Row UL_Row]; 
        ROIColumn2 = [UL_Column UR_Column LR_Column LL_Column UL_Column];
        [Row, Column]= size(A2);
        
        % Region of Interest(ROI) Mask 
        mask= poly2mask(ROIColumn2, ROIRow2, Row, Column);
        ToA2ROI = double(ToA_rad2).*mask;
        
        %Mean of the ToA Radiance of Interested Region
        ToA2Mean(date) = mean(ToA2ROI(ToA2ROI~=0));
        
        %ToA Radiance with a correction for Sun Angle
        Sun_Azimuth(date)= MTL_List.GROUP_IMAGE_ATTRIBUTES.SUN_AZIMUTH;
        Sun_Elevation(date)=MTL_List.GROUP_IMAGE_ATTRIBUTES.SUN_ELEVATION;
        ToA2Mean_corr(date)=ToA2Mean(date)./sind(Sun_Elevation(date));
        
        %-------------------------Band3: Green----------------------------%
        rmb3(date)= MTL_List.RADIOMETRIC_RESCALING.RADIANCE_MULT_BAND_3;
        rab3(date)= MTL_List.RADIOMETRIC_RESCALING.RADIANCE_ADD_BAND_3;
        [A3, R3] = geotiffread(fullfile(base,dates(date).name,'LC1',Band3.name));
        A3=double(A3);
        %DN1(date)=mean2(A1(r,c));
        ToA_rad3=A3*rmb3(date)+rab3(date);  % ToA Radiance calculation                       
      
        % Calculation of pixel coordinates from map coordinates
        [UL_Row, UL_Column] = map2pix(R3, UL_LongitudeL8, UL_LatitudeL8);
        [UR_Row, UR_Column] = map2pix(R3, UR_LongitudeL8, UR_LatitudeL8);
        [LL_Row, LL_Column] = map2pix(R3, LL_LongitudeL8, LL_LatitudeL8);
        [LR_Row, LR_Column] = map2pix(R3, LR_LongitudeL8, LR_LatitudeL8);
        
        ROIRow3 = [UL_Row UR_Row LR_Row LL_Row UL_Row]; 
        ROIColumn3 = [UL_Column UR_Column LR_Column LL_Column UL_Column];
        [Row, Column]= size(A3);
        
        % Region of Interest(ROI) Mask 
        mask= poly2mask(ROIColumn3, ROIRow3, Row, Column);
        ToA3ROI = double(ToA_rad3).*mask;
        
        %Mean of the ToA Radiance of Interested Region
        ToA3Mean(date) = mean(ToA3ROI(ToA3ROI~=0));
        
        %ToA Radiance with a correction for Sun Angle
        Sun_Azimuth(date)= MTL_List.GROUP_IMAGE_ATTRIBUTES.SUN_AZIMUTH;
        Sun_Elevation(date)=MTL_List.GROUP_IMAGE_ATTRIBUTES.SUN_ELEVATION;
        ToA3Mean_corr(date)=ToA3Mean(date)./sind(Sun_Elevation(date));

        %----------------------------Band4:Red----------------------------%
        rmb4(date)= MTL_List.RADIOMETRIC_RESCALING.RADIANCE_MULT_BAND_4;
        rab4(date)= MTL_List.RADIOMETRIC_RESCALING.RADIANCE_ADD_BAND_4;
        [A4, R4] = geotiffread(fullfile(base,dates(date).name,'LC1',Band4.name));
        A4=double(A4);
        %DN1(date)=mean2(A1(r,c));
        ToA_rad4=A4*rmb4(date)+rab4(date);  % ToA Radiance calculation                       
      
        % Calculation of pixel coordinates from map coordinates
        [UL_Row, UL_Column] = map2pix(R4, UL_LongitudeL8, UL_LatitudeL8);
        [UR_Row, UR_Column] = map2pix(R4, UR_LongitudeL8, UR_LatitudeL8);
        [LL_Row, LL_Column] = map2pix(R4, LL_LongitudeL8, LL_LatitudeL8);
        [LR_Row, LR_Column] = map2pix(R4, LR_LongitudeL8, LR_LatitudeL8);
        
        ROIRow4 = [UL_Row UR_Row LR_Row LL_Row UL_Row]; 
        ROIColumn4 = [UL_Column UR_Column LR_Column LL_Column UL_Column];
        [Row, Column]= size(A4);
        
        % Region of Interest(ROI) Mask 
        mask= poly2mask(ROIColumn4, ROIRow4, Row, Column);
        ToA4ROI = double(ToA_rad4).*mask;
        
        %Mean of the ToA Radiance of Interested Region
        ToA4Mean(date) = mean(ToA4ROI(ToA4ROI~=0));
        
        %ToA Radiance with a correction for Sun Angle
        Sun_Azimuth(date)= MTL_List.GROUP_IMAGE_ATTRIBUTES.SUN_AZIMUTH;
        Sun_Elevation(date)=MTL_List.GROUP_IMAGE_ATTRIBUTES.SUN_ELEVATION;
        ToA4Mean_corr(date)=ToA4Mean(date)./sind(Sun_Elevation(date));

        %-----------------Band5: Near Infrared or NIR---------------------%
        rmb5(date)= MTL_List.RADIOMETRIC_RESCALING.RADIANCE_MULT_BAND_5;
        rab5(date)= MTL_List.RADIOMETRIC_RESCALING.RADIANCE_ADD_BAND_5;
        [A5, R5] = geotiffread(fullfile(base,dates(date).name,'LC1',Band5.name));
        A5=double(A5);
        %DN1(date)=mean2(A1(r,c));
        ToA_rad5=A5*rmb5(date)+rab5(date);  % ToA Radiance calculation                       
      
        % Calculation of pixel coordinates from map coordinates
        [UL_Row, UL_Column] = map2pix(R5, UL_LongitudeL8, UL_LatitudeL8);
        [UR_Row, UR_Column] = map2pix(R5, UR_LongitudeL8, UR_LatitudeL8);
        [LL_Row, LL_Column] = map2pix(R5, LL_LongitudeL8, LL_LatitudeL8);
        [LR_Row, LR_Column] = map2pix(R5, LR_LongitudeL8, LR_LatitudeL8);
        
        ROIRow5 = [UL_Row UR_Row LR_Row LL_Row UL_Row]; 
        ROIColumn5 = [UL_Column UR_Column LR_Column LL_Column UL_Column];
        [Row, Column]= size(A5);
        
        % Region of Interest(ROI) Mask 
        mask= poly2mask(ROIColumn5, ROIRow5, Row, Column);
        ToA5ROI = double(ToA_rad5).*mask;
        
        %Mean of the ToA Radiance of Interested Region
        ToA5Mean(date) = mean(ToA5ROI(ToA5ROI~=0));
        
        %ToA Radiance with a correction for Sun Angle
        Sun_Azimuth(date)= MTL_List.GROUP_IMAGE_ATTRIBUTES.SUN_AZIMUTH;
        Sun_Elevation(date)=MTL_List.GROUP_IMAGE_ATTRIBUTES.SUN_ELEVATION;
        ToA5Mean_corr(date)=ToA5Mean(date)./sind(Sun_Elevation(date));
        
        %--------------Band6: Shortwave Infrared1 or SwIR1----------------%
        rmb6(date)= MTL_List.RADIOMETRIC_RESCALING.RADIANCE_MULT_BAND_6;
        rab6(date)= MTL_List.RADIOMETRIC_RESCALING.RADIANCE_ADD_BAND_6;
        [A6, R6] = geotiffread(fullfile(base,dates(date).name,'LC1',Band6.name));
        A6=double(A6);
        %DN1(date)=mean2(A1(r,c));
        ToA_rad6=A6*rmb6(date)+rab6(date);  % ToA Radiance calculation                       
      
        % Calculation of pixel coordinates from map coordinates
        [UL_Row, UL_Column] = map2pix(R6, UL_LongitudeL8, UL_LatitudeL8);
        [UR_Row, UR_Column] = map2pix(R6, UR_LongitudeL8, UR_LatitudeL8);
        [LL_Row, LL_Column] = map2pix(R6, LL_LongitudeL8, LL_LatitudeL8);
        [LR_Row, LR_Column] = map2pix(R6, LR_LongitudeL8, LR_LatitudeL8);
        
        ROIRow6 = [UL_Row UR_Row LR_Row LL_Row UL_Row]; 
        ROIColumn6 = [UL_Column UR_Column LR_Column LL_Column UL_Column];
        [Row, Column]= size(A6);
        
        % Region of Interest(ROI) Mask 
        mask= poly2mask(ROIColumn6, ROIRow6, Row, Column);
        ToA6ROI = double(ToA_rad6).*mask;
        
        %Mean of the ToA Radiance of Interested Region
        ToA6Mean(date) = mean(ToA6ROI(ToA6ROI~=0));
        
        %ToA Radiance with a correction for Sun Angle
        Sun_Azimuth(date)= MTL_List.GROUP_IMAGE_ATTRIBUTES.SUN_AZIMUTH;
        Sun_Elevation(date)=MTL_List.GROUP_IMAGE_ATTRIBUTES.SUN_ELEVATION;
        ToA6Mean_corr(date)=ToA6Mean(date)./sind(Sun_Elevation(date));

        %----------------Band7: Shortwave Infrared or SwIR2---------------%
        rmb7(date)= MTL_List.RADIOMETRIC_RESCALING.RADIANCE_MULT_BAND_7;
        rab7(date)= MTL_List.RADIOMETRIC_RESCALING.RADIANCE_ADD_BAND_7;
        [A7, R7] = geotiffread(fullfile(base,dates(date).name,'LC1',Band7.name));
        A7=double(A7);
        %DN7(date)=mean2(A1(r,c));
        ToA_rad7=A7*rmb7(date)+rab7(date);  % ToA Radiance calculation                       
      
        % Calculation of pixel coordinates from map coordinates
        [UL_Row, UL_Column] = map2pix(R7, UL_LongitudeL8, UL_LatitudeL8);
        [UR_Row, UR_Column] = map2pix(R7, UR_LongitudeL8, UR_LatitudeL8);
        [LL_Row, LL_Column] = map2pix(R7, LL_LongitudeL8, LL_LatitudeL8);
        [LR_Row, LR_Column] = map2pix(R7, LR_LongitudeL8, LR_LatitudeL8);
        
        ROIRow7 = [UL_Row UR_Row LR_Row LL_Row UL_Row]; 
        ROIColumn7 = [UL_Column UR_Column LR_Column LL_Column UL_Column];
        [Row, Column]= size(A7);
        
        % Region of Interest(ROI) Mask 
        mask= poly2mask(ROIColumn7, ROIRow7, Row, Column);
        ToA7ROI = double(ToA_rad7).*mask;
        
        %Mean of the ToA Radiance of Interested Region
        ToA7Mean(date) = mean(ToA7ROI(ToA7ROI~=0));
        
        %ToA Radiance with a correction for Sun Angle
        Sun_Azimuth(date)= MTL_List.GROUP_IMAGE_ATTRIBUTES.SUN_AZIMUTH;
        Sun_Elevation(date)=MTL_List.GROUP_IMAGE_ATTRIBUTES.SUN_ELEVATION;
        ToA7Mean_corr(date)=ToA7Mean(date)./sind(Sun_Elevation(date));
        
        
      end

   
  end
  
  % Removing Outliers
  
% s1= size(ToA1Mean_corr);
% s2= size(ToA2Mean_corr);
% s3= size(ToA3Mean_corr);
% s4= size(ToA4Mean_corr);
% s5= size(ToA5Mean_corr);
% s6= size(ToA6Mean_corr);
% s7= size(ToA7Mean_corr);
% 
% outlier1 = abs(ToA1Mean_corr(:,(1:s1(:,2)))-mean(ToA1Mean_corr(:,(1:s1(:,2))))>2*std(ToA1Mean_corr(:,(1:s1(:,2)))));
% outlier_index1 = find(outlier1);

% outlier2 = abs(ToA2Mean_corr(:,(1:s2(:,2)))-mean(ToA2Mean_corr(:,(1:s2(:,2))))>2*std(ToA2Mean_corr(:,(1:s2(:,2)))));
% outlier_index2 = find(outlier2);
% outlier3 = abs(ToA3Mean_corr(:,(1:s3(:,2)))-mean(ToA3Mean_corr(:,(1:s3(:,2))))>2*std(ToA3Mean_corr(:,(1:s3(:,2)))));
% outlier_index3 = find(outlier3);
% outlier4 = abs(ToA4Mean_corr(:,(1:s4(:,2)))-mean(ToA4Mean_corr(:,(1:s4(:,2))))>2*std(ToA4Mean_corr(:,(1:s4(:,2)))));
% outlier_index4 = find(outlier4);
% outlier5 = abs(ToA5Mean_corr(:,(1:s5(:,2)))-mean(ToA5Mean_corr(:,(1:s5(:,2))))>2*std(ToA5Mean_corr(:,(1:s5(:,2)))));
% outlier_index5 = find(outlier5);
% outlier6 = abs(ToA6Mean_corr(:,(1:s6(:,2)))-mean(ToA6Mean_corr(:,(1:s6(:,2))))>2*std(ToA6Mean_corr(:,(1:s6(:,2)))));
% outlier_index6 = find(outlier6);
% outlier7 = abs(ToA7Mean_corr(:,(1:s7(:,2)))-mean(ToA7Mean_corr(:,(1:s7(:,2))))>2*std(ToA7Mean_corr(:,(1:s7(:,2)))));
% outlier_index7 = find(outlier7);
% 
% ToA1Mean_corr(outlier_index1) = [];
% ToA2Mean_corr(outlier_index2) = [];
% ToA3Mean_corr(outlier_index3) = [];
% ToA4Mean_corr(outlier_index4) = [];
% ToA5Mean_corr(outlier_index5) = [];
% ToA6Mean_corr(outlier_index6) = [];
% ToA7Mean_corr(outlier_index7) = [];
% 
% DoY(outlier_index1)=[];
% DeciYear(outlier_index1)=[];
