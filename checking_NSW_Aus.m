% clc, close, clear all;
clearvars -except all_mean all_RSD all_SD

% % Path/row-090/082- New South Wells, Australia
base='Z:\SpecialNeeds\Sakib\DoveData\P090\R082\';
dates = dir(base);dates([1 2])=[];

equator= 10000000;

for ROI=1:4
    
  if ROI==1
    % ROI Coordinates
    UL_x=216608;
    UL_y=6430697-equator;
    UR_x=217096;
    UR_y=6430822-equator;
    LR_x=217380;
    LR_y=6430116-equator;
    LL_x=216863;
    LL_y=6429910-equator;

    
        elseif ROI==2
    UL_x=230132;
    UL_y=6432712-equator;
    UR_x=230542;
    UR_y=6432656-equator;
    LR_x=230525;
    LR_y=6432476-equator;
    LL_x=230124;
    LL_y=6432558-equator;

       elseif ROI==3
    UL_x=228172;
    UL_y=6433448-equator;
    UR_x=228342;
    UR_y=6433428-equator;
    LR_x=228330;
    LR_y=6433212-equator;
    LL_x=228144;
    LL_y=6433244-equator;



      elseif ROI==4

        % ROI Coordinates
    UL_x=228360;
    UL_y=6431255-equator;
    UR_x=228464;
    UR_y=6431203-equator;
    LR_x=228441;
    LR_y=6431134-equator;
    LL_x=228333;
    LL_y=6431188-equator;


%       elseif ROI==5
%     UL_x=671265;
%     UL_y=6405755-equator;
%     UR_x=671517;
%     UR_y=6405669-equator;
%     LR_x=671416;
%     LR_y=6405591-equator;
%     LL_x=671229;
%     LL_y=6405666-equator;

  end
    
x_vec=[UL_x UR_x LR_x LL_x UL_x];
y_vec=[UL_y UR_y LR_y LL_y UL_y];

bands=  {'1' '2' '3' '4' '5' '6' '7'}; 
band_str= {'CA','B',  'G', 'R', 'NIR', 'SWIR1', 'SWIR2'};
band_name = {'CA' ,'BLUE','GREEN' ,'RED' ,'NIR' ,'SWIR1' ,'SWIR2'};
band_colors={'b','c','g','r','m','[0.6667 0 0]','k'};

%MTL File Parsing 
MTL=dir(fullfile(base, dates.name,'*MTL.txt'));
[MTL_List_L8, vaule]= MTL_parser_L8(fullfile(MTL.folder, MTL.name));

%Extracting all the values from Metadata file
MData_file=strcat(base, filesep, dates.name, filesep, '20181111_002927_91_1047_3B_AnalyticMS_metadata.xml');
[MData_values]= xml2struct_new_v(MData_file);
[MData_values_old]= xml2struct(MData_file);
%%
for D_band=1:4 %dove band
    %Image file location
    L8_band=D_band+1;
    Folder_info=dir(fullfile(base, dates.name, strcat('*B', bands{L8_band},'.tif')));  
    L8_image_file= fullfile(Folder_info.folder, Folder_info.name);
    %Dove_image_file= strcat(base,filesep, dates.name, filesep, '002927_91_cubicspline.tif');
    Dove_image_file= strcat(base,filesep, dates.name, filesep, '20181111_002927_91_1047_3B_AnalyticMS.tif');
     
    %Reading the image
    [DN, R_L8] = geotiffread(L8_image_file);
    DN=double(DN);
   
    L8_info = geotiffinfo(L8_image_file);
    %[Dove_image_all_band, R_dove]=geotiffread(Dove_image_file);
    lat_TL = str2num(MData_values.ps_colon_EarthObservation.gml_colon_target...
            .ps_colon_Footprint.ps_colon_geographicLocation.ps_colon_topLeft.ps_colon_latitude.Text);   
    lon_TL = str2num(MData_values.ps_colon_EarthObservation.gml_colon_target...
            .ps_colon_Footprint.ps_colon_geographicLocation.ps_colon_topLeft.ps_colon_longitude.Text);
     
    [X, Y] = ll2utm(lat_TL, lon_TL);
    R_dove = [0 -3; 3 0; 205194 Y-equator];
    Dove_image_all_band =imread(Dove_image_file);
    
      %%  
    %%%% For pixel wise SNR
     close all
     band_name = {'BLUE','GREEN' ,'RED' ,'NIR' };
     band_colors={'c','g','r','m'};
    
    for D_band=1:4
        % all the Angles
        Dove_image_info.D_SEle_Angle=str2num(MData_values.ps_colon_EarthObservation.gml_colon_using...
                .eop_colon_EarthObservationEquipment.eop_colon_acquisitionParameters.ps_colon_Acquisition...
                .opt_colon_illuminationElevationAngle.Text);  
        Dove_image_info.D_SAzi_Angle=str2num(MData_values.ps_colon_EarthObservation.gml_colon_using...
                .eop_colon_EarthObservationEquipment.eop_colon_acquisitionParameters.ps_colon_Acquisition...
                .opt_colon_illuminationAzimuthAngle.Text); 

        Dove_image=Dove_image_all_band(:,:,D_band); %Blue band
        Dove_image=double(Dove_image);
        DN_Dove=double(Dove_image);
        DN_Dove(DN_Dove==0)=nan;

        band_sf=0.01;
        TOA_rad_D = DN_Dove.*band_sf;
        TOA_rad_D = TOA_rad_D./sind(Dove_image_info.D_SEle_Angle);

        i = fspecial('average', [91 91]);
        Dove_mean_filter = filter2(i, TOA_rad_D);
        Dove_SD_filter = stdfilt(TOA_rad_D);
        SNR_esti = Dove_mean_filter./Dove_SD_filter;
        
        MeanTOA_rad_D=(Dove_mean_filter(~isnan(Dove_mean_filter)));
        SNR_esti_D=(SNR_esti(~isnan(SNR_esti))); 
        plot(MeanTOA_rad_D(1:10:end), SNR_esti_D(1:10:end), '.')
        
%         x = decimate(MeanTOA_rad_D,4);
%         y = decimate(SNR_esti_D,4);
        subplot(2,2,D_band)
        plot(MeanTOA_rad_D, SNR_esti_D , 'color', band_colors{D_band}, 'LineStyle','None','Marker','.','markers', 10)
        %ylim([0 400])
        hold on
        
        %%% Dove   
        % Radiance
        title(strcat(strcat(band_name{D_band}, ' Band')));
        xlabel('mean TOA Radiance of 1047')
        ylabel('SNR')
        
        hold on
        grid on
        grid minor
        ax  = gca;
        ax.FontSize = 8;
        ax.GridColor = 'k';
    %   ax.GridAlpha = 0.8;
        ax.MinorGridColor = 'k';
        %ax.FontName = 'Times New Roman';
    end
    

    
    
    % Mulitplicative and additive factors
    rmb= MTL_List_L8.RADIOMETRIC_RESCALING.(strcat('RADIANCE_MULT_BAND_', bands{L8_band}));
    rab= MTL_List_L8.RADIOMETRIC_RESCALING.(strcat('RADIANCE_ADD_BAND_', bands{L8_band}));
    
    % Mulitplicative and additive factors
    refmb= MTL_List_L8.RADIOMETRIC_RESCALING.(strcat('REFLECTANCE_MULT_BAND_', bands{L8_band}));
    refab= MTL_List_L8.RADIOMETRIC_RESCALING.(strcat('REFLECTANCE_ADD_BAND_', bands{L8_band}));
    
    %Map coordinate to pixel coordinate
    [Pixel_Row_unrounded_L8, Pixel_Column_unrounded_L8] = map2pix(R_L8, x_vec, y_vec);
    Pixel_Row_L8= round(Pixel_Row_unrounded_L8);
    Pixel_Column_L8= round(Pixel_Column_unrounded_L8);
    [Row_L8, Column_L8]= size(DN);

    % ROI Mask 
    mask_L8= poly2mask(Pixel_Column_L8, Pixel_Row_L8,  Row_L8, Column_L8);
    DN=DN.*mask_L8;
    DN(DN==0)= NaN;
    L8_TOArad=DN*rmb+rab;
  
    L8_TOAref=DN*refmb+refab;
    L8_TOArad(L8_TOArad<0)=nan;
    
    L8_TOAref(L8_TOAref<0)=nan;
    L8_TOAref(L8_TOAref==0)=nan;
     
    %%%%%%%%%%%%Dove 
    Dove_image=Dove_image_all_band(:,:,D_band); %Blue band
    Dove_image=double(Dove_image);
    DN_Dove=double(Dove_image);
    DN_Dove(DN_Dove==0)=nan;
     
    %Map coordinate to pixel coordinate
    [Pixel_Row_unrounded_d, Pixel_Column_unrounded_d] = map2pix(R_dove, x_vec, y_vec);
    Pixel_Row_d= round(Pixel_Row_unrounded_d);
    Pixel_Column_d= round(Pixel_Column_unrounded_d);

    [Row_d, Column_d]= size(DN_Dove);
    mask_d= poly2mask(Pixel_Column_d, Pixel_Row_d, Row_d, Column_d);
    Image_dove_masked=DN_Dove.*mask_d;

      % all the Angles
    Dove_image_info.D_SEle_Angle=str2num(MData_values.ps_colon_EarthObservation.gml_colon_using...
    .eop_colon_EarthObservationEquipment.eop_colon_acquisitionParameters.ps_colon_Acquisition...
    .opt_colon_illuminationElevationAngle.Text);  
    Dove_image_info.D_SAzi_Angle=str2num(MData_values.ps_colon_EarthObservation.gml_colon_using...
    .eop_colon_EarthObservationEquipment.eop_colon_acquisitionParameters.ps_colon_Acquisition...
    .opt_colon_illuminationAzimuthAngle.Text);  
   
    band_sf=0.01;
    D_ROIrad=Image_dove_masked.*(band_sf);
    D_ROIrad(D_ROIrad==0)=nan;
    D_ROIrad=D_ROIrad./sind(Dove_image_info.D_SEle_Angle);
    
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
    Ref_coeff =str2num(MData_values_old.Children(10).Children(2).Children(id).Children(10).Children.Data);
     
    D_ROIref=Image_dove_masked.*(Ref_coeff);
   % D_ROIref=Image_dove_masked.*(Ref_coeff);
    D_ROIref(D_ROIref==0)=nan;
    D_ROIref(D_ROIref<0)=nan;
     
   % D_ROIref=D_ROIref./sind(Dove_image_info.D_SEle_Angle);
    
   % DN_D_Mean=nanmean(nanmean(DN_Dove));
    
    %D_ROIrad_sampled=imresize(D_ROIrad, 3/30); % sampling
    
    % Solar angles- scene center
    Sun_Azimuth_L8= MTL_List_L8.GROUP_IMAGE_ATTRIBUTES.SUN_AZIMUTH;
    Sun_Elevation_L8=MTL_List_L8.GROUP_IMAGE_ATTRIBUTES.SUN_ELEVATION;
   
    %L8_ROIrad_b2=L8_ROIrad_b2./sind(Sun_Elevation); %cosine correction

    % New corner Coordinates
%     [nR,nC]=map2pix(R_L8, R_dove.XWorldLimits(1,1), R_dove.YWorldLimits(1,2));
%     nR=round(nR);
%     nC=round(nC);
%     [sR, sC]= size(D_ROIrad);
%     D_TOArad_final = NaN(size(DN,1), size(DN,2));%size is same for all L8 band
%     D_TOArad_final(nR:nR+sR-1,nC:nC+sC-1)= D_ROIrad;
%     
%     D_TOAref_final = NaN(size(DN,1), size(DN,2));%size is same for all L8 band
%     D_TOAref_final(nR:nR+sR-1,nC:nC+sC-1)= D_ROIref;
%     D_TOAref_final(D_TOAref_final<0)= NaN;
%     D_TOAref_final(D_TOAref_final== 0)= NaN;
    
    D_TOArad_final = D_ROIrad;
    D_TOAref_final = D_ROIref;
    
    %L8_TOArad_final = L8_TOArad;
    
    % Angles for cosine correction
    Angle_folder= dir(fullfile(base, dates.name, strcat('*solar_B05.img')));
    L8_Angle_file= fullfile( Angle_folder.folder, Angle_folder.name);
    SolarAngleInfo = multibandread(L8_Angle_file,[size(DN),2],'int16',0,'bsq','ieee-le');
    %solar_azimuth=(SolarAngleInfo(:,:,1))/100;
    solar_zenith=(SolarAngleInfo(:,:,2))/100;
    
    %L8 solar zenith pixel by pixel
    L8_mat_logical = ~isnan(L8_TOArad);
    solar_zenith_L8ROImat = solar_zenith.*L8_mat_logical;
    solar_zenith_L8ROImat(solar_zenith_L8ROImat==0)=nan;
     
    L8_TOArad_final = L8_TOArad./cosd(solar_zenith_L8ROImat); % check is it in degree or radian 
    L8_TOArad_final( L8_TOArad_final<0)= NaN;
    L8_TOArad_final( L8_TOArad_final== 0)= NaN;
    
    L8_TOAref_final = L8_TOAref./cosd(solar_zenith_L8ROImat);
    L8_TOAref_final(L8_TOAref_final<0)= NaN;
    L8_TOAref_final(L8_TOAref_final == 0)= NaN;
    
    %using scene center angle
    L8_TOArad_final2 = L8_TOArad./sind(Sun_Elevation_L8);
    L8_TOAref_final2 = L8_TOAref./sind(Sun_Elevation_L8);
    
     %Storing radiances to different varibles
%     if L8_band==2
%         L8_TOArad_band2=L8_TOArad_final;
%         D_TOArad_band1=D_TOArad_final;
%         
%         Mean_TOArad_DROI_b1 = nanmean(nanmean(D_TOArad_band1));
%         Mean_TOArad_L8ROI_b2 = nanmean(nanmean(L8_TOArad_band2));
%        % clearvars L8_TOArad_final D_TOArad_final
%         
%     elseif L8_band==3
%         L8_TOArad_band3=L8_TOArad_final;
%         D_TOArad_band2=D_TOArad_final;
%         
%         Mean_TOArad_DROI_b2 = nanmean(nanmean(D_TOArad_band2));
%         Mean_TOArad_L8ROI_b3 = nanmean(nanmean(L8_TOArad_band3));
%        % clearvars L8_TOArad_final D_TOArad_final
%         
%     elseif L8_band==4
%         L8_TOArad_band4=L8_TOArad_final;
%         D_TOArad_band3=D_TOArad_final;
%      
%         Mean_TOArad_DROI_b3 = nanmean(nanmean(D_TOArad_band3));
%         Mean_TOArad_L8ROI_b4 = nanmean(nanmean(L8_TOArad_band4));
%         %clearvars L8_TOArad_final D_TOArad_final
%         
%     elseif L8_band==5
%         L8_TOArad_band5=L8_TOArad_final;
%         D_TOArad_band4=D_TOArad_final;
%         Mean_TOArad_DROI_b4 = nanmean(nanmean(D_TOArad_band4));
%         Mean_TOArad_L8ROI_b5 = nanmean(nanmean(L8_TOArad_band5));
%         %clearvars L8_TOArad_final D_TOArad_final
%     end
        %%% Radiance
        temp_L8_rad = L8_TOArad_final;
        temp_L8_rad=(temp_L8_rad(~isnan(temp_L8_rad)));

        temp_D_rad = D_TOArad_final;
        temp_D_rad=(temp_D_rad(~isnan(temp_D_rad)));

        % Mean
        mean_D_TOArad(D_band)= mean(temp_D_rad);
        mean_L8_TOArad(D_band)= mean(temp_L8_rad); 

        % Standard Deviation
        Std_D_TOArad(D_band)= std(temp_D_rad);
        Std_L8_TOArad(D_band)= std(temp_L8_rad);

        %%% Reflectance
        temp_D_ref = D_TOAref_final;
        temp_D_ref=(temp_D_ref(~isnan(temp_D_ref)));

        temp_L8_ref = L8_TOAref_final;
        temp_L8_ref=(temp_L8_ref(~isnan(temp_L8_ref)));

        % Mean
        mean_D_TOAref(D_band)= mean(temp_D_ref);
        mean_L8_TOAref(D_band)= mean(temp_L8_ref);

        %Standard Deviation
        Std_D_TOAref(D_band)= std(temp_D_ref);
        Std_L8_TOAref(D_band)= std(temp_L8_ref);
end
  
%     all_mean.D.ROI(ROI,:) = mean_D_TOArad;
%     all_mean.L8.ROI(ROI,:) = mean_L8_TOArad;
    
     all_mean_NSW_Aus2.D.ROI(ROI,:) = mean_D_TOArad;
     all_mean_NSW_Aus2.L8.ROI(ROI,:) = mean_L8_TOArad;
    % all_mean.L8_2.ROI(ROI,:)= mean_L8_2_TOArad;
   
     all_mean_ref_NSW_Aus2.D.ROI(ROI,:) = mean_D_TOAref;
     all_mean_ref_NSW_Aus2.L8.ROI(ROI,:) = mean_L8_TOAref;
     
      %Standard Deviation- Radiance
     all_SD_NSW_Aus2.D.ROI(ROI,:)= Std_D_TOArad;
     all_SD_NSW_Aus2.L8.ROI(ROI,:)= Std_L8_TOArad;
     
     %Standard Deviation- Reflectance
     all_SD_ref_NSW_Aus2.D.ROI(ROI,:)= Std_D_TOAref;
     all_SD_ref_NSW_Aus2.L8.ROI(ROI,:)= Std_L8_TOAref;
     
     %Relative Standard Deviation from Radiance
     all_RSD_NSW_Aus2.D.ROI(ROI,:) = Std_D_TOArad./mean_D_TOArad;
     all_RSD_NSW_Aus2.L8.ROI(ROI,:) = Std_L8_TOArad./mean_L8_TOArad;
     
     %Relative Standard Deviation from Reflectance
     all_RSD_ref_NSW_Aus2.D.ROI(ROI,:) = Std_D_TOAref./mean_D_TOAref;
     all_RSD_ref_NSW_Aus2.L8.ROI(ROI,:) = Std_L8_TOAref./mean_L8_TOAref; 
end

    all_mean.Radiance.n_NSW_Aus2 = all_mean_NSW_Aus2;
    all_mean.Reflectance.n_NSW_Aus2 = all_mean_ref_NSW_Aus2;
    
    all_SD.Radiance.n_NSW_Aus2 = all_SD_NSW_Aus2;
    all_SD.Reflectance.n_NSW_Aus2 = all_SD_ref_NSW_Aus2;
    
    all_RSD.Radiance.n_NSW_Aus2 = all_RSD_NSW_Aus2;
    all_RSD.Reflectance.n_NSW_Aus2= all_RSD_ref_NSW_Aus2;

%%
%%%%%%%% Pixel-by-Pixel Plot
%Pixel-by-Pixel Comparison of TOA Radiance of L8 and Dove
close all
band_name = {'BLUE','GREEN' ,'RED' ,'NIR' };
band_colors={'c','g','r','m'};
leg ={'ROI 1', 'ROI 2', 'ROI 3',  'ROI 4', 'ROI 5',  'ROI 6', 'ROI 7'};
marker = {'o','s', '^', 'h', '>','p', 'd'};
 for band=1:4 %Dove band
     
     if band==1
         subplot(2,2,1); 
%         plot(all_mean.L8.ROI(:,band), all_mean.D.ROI(:,band), 'color', band_colors{band}, 'LineStyle','None','Marker','.','markers', 18)
%         hold on
        L8_band2 = all_mean.L8.ROI(:,band);
        D_band1 = all_mean.D.ROI(:,band);
       % figure(1)
        hold on
        sz = 120;
        for k=1:length(all_mean.L8.ROI(:,band))
             scatter(L8_band2(k), D_band1(k),sz, marker{k}, 'filled', band_colors{band})
        end
        hold off
        legend(leg,'FontSize',16, 'Location','southeast');
        
        hold on
        plot([0 100], [0 100], 'k')
        
     elseif band==2
        subplot(2,2,2); 
%         plot(all_mean.L8.ROI(:,band), all_mean.D.ROI(:,band), 'color', band_colors{band}, 'LineStyle','None','Marker','.','markers', 18)
        L8_band3 = all_mean.L8.ROI(:,band);
        D_band2 = all_mean.D.ROI(:,band);
        
       % figure(2)
        hold on
        sz = 120;
        for k=1:length(all_mean.L8.ROI(:,band))
             scatter(L8_band3(k), D_band2(k),sz, marker{k}, 'filled', band_colors{band})
        end
        hold off
        legend(leg,'FontSize',16, 'Location','southeast');
        
        hold on
        plot([0 100], [0 100], 'k')
        
     elseif band==3
       subplot(2,2,3);
%         plot(all_mean.L8.ROI(:,band), all_mean.D.ROI(:,band), 'color', band_colors{band}, 'LineStyle','None','Marker','.','markers', 18)

        L8_band4 = all_mean.L8.ROI(:,band);
        D_band3 = all_mean.D.ROI(:,band);
      %  figure(3)
        hold on
        sz = 120;
        
        for k=1:length(all_mean.L8.ROI(:,band))
             scatter(L8_band4(k), D_band3(k),sz, marker{k}, 'filled', band_colors{band})
        end
        
        hold off
        
        legend(leg,'FontSize',16, 'Location','southeast');
        hold on
         plot([0 100], [0 100], 'k')
        
     elseif band==4
        subplot(2,2,4);
%         plot(all_mean.L8.ROI(:,band), all_mean.D.ROI(:,band), 'color', band_colors{band}, 'LineStyle','None','Marker','.','markers', 18)
          
        L8_band5 = all_mean.L8.ROI(:,band);
        D_band4 = all_mean.D.ROI(:,band);
        
        %figure(4)
        hold on
        sz = 120;
        
        for k=1:length(all_mean.L8.ROI(:,band))
             scatter(L8_band5(k), D_band4(k),sz, marker{k}, 'filled', band_colors{band})
        end
        
        hold off
        legend(leg,'FontSize',16, 'Location','southeast');
        
        hold on
        plot([0 100], [0 100], 'k')
        
     end
     
    title(strcat('Mean TOA Radiance Comparison of L8 and Dove', {', '}, strcat(band_name{band}, ' Band')));
    xlabel('Mean TOA Radiance L8 ROI (W/Sr/m^2/{\mum})')
    ylabel('Mean TOA Radiance 1047 ROI (W/Sr/m^2/{\mum})')

    hold on
    grid on
    grid minor
    ax  = gca;
    ax.FontSize = 12;
    ax.GridColor = 'k';
%   ax.GridAlpha = 0.8;
    ax.MinorGridColor = 'k';
    %ax.FontName = 'Times New Roman';
 end
%%
 
 
% for ROI=1:2
%     all_mean.D.ROI(ROI,:) = mean_D_TOArad;
%     all_mean.D.ROI(ROI,:) = mean_D_TOArad;
%     
% end
% % base
% BQA=geotiffread('Z:\ImageDrive\OLI.TIRS\L8\P180\R040\20150314\LC1\LC08_L1TP_180040_20150314_20170412_01_T1_BQA.tif');
 SolarAngleInformation = multibandread('Z:\ImageDrive\OLI.TIRS\L8\P180\R040\20150314\LC1\LC08_L1TP_180040_20150314_20170412_01_T1_solar_B05.img',[size(DN),2],'int16',0,'bsq','ieee-le');
 solar_zenith=(SolarAngleInformation(:,:,2))/100;
 solar_azimuth=(SolarAngleInformation(:,:,1))/100;

% %%
% all_mean.D=Mean_TOArad_DROI_b3;
% all_mean.L8=Mean_TOArad_L8ROI_b4;
% all_mean_temp=all_mean;