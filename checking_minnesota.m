clc, close, clear all;
%clearvars -except all_mean;

% % Path/row-029/028- Minnesota
base='Z:\SpecialNeeds\Sakib\DoveData\P029\R028\';

dates = dir(base);
dates([1 2])=[];

equator= 10000000;

for ROI=1:4
    
  if ROI==1
    % ROI Coordinates
    UL_x=272459;
    UL_y=5016702;
    UR_x=272992;
    UR_y=5016676;
    LR_x=272968;
    LR_y=5016174;
    LL_x=272436;
    LL_y=5016204;


        elseif ROI==2

    UL_x=277221;
    UL_y=5015720;
    UR_x=277384;
    UR_y=5015717;
    LR_x=277443;
    LR_y=5015169;
    LL_x=277198;
    LL_y=5015161;


       elseif ROI==3
    UL_x=277547;
    UL_y=5023661;
    UR_x=278116;
    UR_y=5023632;
    LR_x=278047;
    LR_y=5023132;
    LL_x=277547;
    LL_y=5023132;


      elseif ROI==4

        % ROI Coordinates
    UL_x=282508;
    UL_y=5026826;
    UR_x=283093;
    UR_y=5026792;
    LR_x=283045;
    LR_y=5026241;
    LL_x=282486;
    LL_y=5026278;


%     elseif ROI==5
% 
%         % ROI Coordinates
%     UL_x=656294;
%     UL_y=5004222;
%     UR_x=656486;
%     UR_y=5004230;
%     LR_x=656508;
%     LR_y=5003843;
%     LL_x=656341;
%     LL_y=5003847;
    

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
MData_file=strcat(base, dates.name, filesep, '20181217_160118_89_1047_3B_AnalyticMS_metadata.xml');
[MData_values]= xml2struct_new_v(MData_file);

for D_band=1:4 %dove band
    %Image file location
    L8_band=D_band+1;
    Folder_info=dir(fullfile(base, dates.name, strcat('*B', bands{L8_band},'.tif')));  
    L8_image_file= fullfile(Folder_info.folder, Folder_info.name);
    Dove_image_file= strcat(base,filesep, dates.name, filesep, '160118_89_cubicspline.tif');
    
    %Reading the image
    [DN, R_L8] = geotiffread(L8_image_file);
    DN=double(DN);
    L8_info = geotiffinfo(L8_image_file);
    [Dove_image_all_band, R_dove]=geotiffread(Dove_image_file);
    
    % Mulitplicative and additive factors
    rmb= MTL_List_L8.RADIOMETRIC_RESCALING.(strcat('RADIANCE_MULT_BAND_', bands{L8_band}));
    rab= MTL_List_L8.RADIOMETRIC_RESCALING.(strcat('RADIANCE_ADD_BAND_', bands{L8_band}));
    
    %Map coordinate to pixel coordinate
    [Pixel_Row_unrounded_L8, Pixel_Column_unrounded_L8] = map2pix(R_L8, x_vec, y_vec);
    Pixel_Row_L8= round(Pixel_Row_unrounded_L8);
    Pixel_Column_L8= round(Pixel_Column_unrounded_L8);
    [Row_L8, Column_L8]= size(DN);

    % ROI Mask 
    mask_L8= poly2mask(Pixel_Column_L8, Pixel_Row_L8,  Row_L8, Column_L8);
    DN=DN.*mask_L8;
    L8_TOArad=DN*rmb+rab;
    
    %%%%%%%%%%%%Dove 
    Dove_image=Dove_image_all_band(:,:,D_band); %Blue band
    Dove_image=double(Dove_image);
    
    %Map coordinate to pixel coordinate
    [Pixel_Row_unrounded_d, Pixel_Column_unrounded_d] = map2pix(R_dove, x_vec, y_vec);
    Pixel_Row_d= round(Pixel_Row_unrounded_d);
    Pixel_Column_d= round(Pixel_Column_unrounded_d);

    [Row_d, Column_d]= size(Dove_image);
    mask_d= poly2mask(Pixel_Column_d, Pixel_Row_d, Row_d, Column_d);
    Image_dove_masked=Dove_image.*mask_d;

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
    
    %D_ROIrad_sampled=imresize(D_ROIrad, 3/30); % sampling
    
    % Solar angles- scene center
    Sun_Azimuth_L8= MTL_List_L8.GROUP_IMAGE_ATTRIBUTES.SUN_AZIMUTH;
    Sun_Elevation_L8=MTL_List_L8.GROUP_IMAGE_ATTRIBUTES.SUN_ELEVATION;
   
    %L8_ROIrad_b2=L8_ROIrad_b2./sind(Sun_Elevation); %cosine correction
    L8_TOArad(L8_TOArad<0)=nan;

    % New corner Coordinates
    [nR,nC]=map2pix(R_L8, R_dove.XWorldLimits(1,1), R_dove.YWorldLimits(1,2));
    nR=round(nR);
    nC=round(nC);
    [sR, sC]= size(D_ROIrad);
    D_TOArad_final = NaN(size(DN,1), size(DN,2));%size is same for all L8 band
    D_TOArad_final(nR:nR+sR-1,nC:nC+sC-1)= D_ROIrad;
    
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
    L8_TOArad_final2 = L8_TOArad./sind(36.51);
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
    mean_D_TOArad(D_band)= nanmean(nanmean(D_TOArad_final));
    mean_L8_TOArad(D_band)= nanmean(nanmean(L8_TOArad_final));
    mean_L8_2_TOArad(D_band)= nanmean(nanmean(L8_TOArad_final2));
    % clearvars L8_TOArad_final D_TOArad_final
end
  
%     all_mean.D.ROI(ROI,:) = mean_D_TOArad;
%     all_mean.L8.ROI(ROI,:) = mean_L8_TOArad;
    
     all_mean.D.ROI(ROI,:) = mean_D_TOArad;
     all_mean.L8.ROI(ROI,:) = mean_L8_TOArad;
     all_mean.L8_2.ROI(ROI,:)= mean_L8_2_TOArad;
   
end

%%
%%%%%%%% Pixel-by-Pixel Plot
%Pixel-by-Pixel Comparison of TOA Radiance of L8 and Dove
close all
band_name = {'BLUE','GREEN' ,'RED' ,'NIR' };
band_colors={'c','g','r','m'};
leg ={'ROI 1', 'ROI 2', 'ROI 3',  'ROI 4', 'ROI 5',  'ROI 6', 'ROI 7'};
marker = {'o', 'p', '^', 'h', '>', 's', 'd'};
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
        plot([0 180], [0 180], 'k')
        
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