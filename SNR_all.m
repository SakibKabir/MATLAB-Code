 %%%% For pixel wise SNR
 clear all
 band_name = {'BLUE', 'GREEN', 'RED', 'NIR'};
 band_colors={'c','g','r','m'};
 bands=  {'1' '2' '3' '4' '5' '6' '7'}; 
 leg_location = {'a_California', 'b_Wisconsin', 'c_Indiana', 'd_Missouri', 'e_J_SA1', 'f_J_SA2',...
                'g_J_SA3', 'h_J_SA4', 'i_CT_SA_1', 'j_CT2_SA_2', 'k_Auckland_NZ', 'l_Sydney_Australia',...
                'm_NSW_Aus1', 'n_NSW_Aus2', 'o_Melbourne_Australia', 'p_Libya4_S52', 'q_Libya4_S50', 'r_illinois',...
                's_Minnesota', 't_SDakota'};
            
for location = 1:20
     
     if location == 1 % California
         base='Z:\SpecialNeeds\Sakib\DoveData\P039\R037\';  
         dates = dir(base); dates([1 2])=[];
         Dove_image_file= strcat(base, dates.name, filesep, '20181223_172655_38_1047_3B_AnalyticMS.tif');
         
         MData_file=strcat(base, filesep, dates.name, filesep, '20181223_172655_38_1047_3B_AnalyticMS_metadata.xml');
         [MData_values]= xml2struct_new_v(MData_file);
         [MData_values_old]= xml2struct(MData_file);
         
     elseif location == 2 % Wisconsin
         base='Z:\SpecialNeeds\Sakib\DoveData\P026\R029\';  
         dates = dir(base); dates([1 2])=[];
         Dove_image_file= strcat(base, dates.name, filesep, '20180923_155547_1047_3B_AnalyticMS.tif');
         
     elseif location == 3 % Indiana
         base='Z:\SpecialNeeds\Sakib\DoveData\P021\R032\';
         dates = dir(base); dates([1 2])=[];  
         Dove_image_file= strcat(base, dates.name, filesep, '20181107_152721_39_1047_3B_AnalyticMS.tif');

     elseif location == 4 % Missouri
         base='Z:\SpecialNeeds\Sakib\DoveData\P023\R034\';
         dates = dir(base); dates([1 2])=[];
         Dove_image_file= strcat(base, dates.name, filesep, '20181121_154801_04_1047_3B_AnalyticMS.tif'); 

     elseif location == 5 % J_SA1
         base='Z:\SpecialNeeds\Sakib\DoveData\P170\R078\';
         dates = dir(base); dates([1 2])=[];
         Dove_image_file= strcat(base, dates.name, filesep, '20181213_082737_49_1047_3B_AnalyticMS.tif');

     elseif location == 6 % J_SA2
         base='Z:\SpecialNeeds\Sakib\DoveData\P170\R079\';
         dates = dir(base); dates([1 2])=[];
         Dove_image_file= strcat(base, dates.name, filesep, '20181213_082717_69_1047_3B_AnalyticMS.tif');

     elseif location == 7 % J_SA3
         base='Z:\SpecialNeeds\Sakib\DoveData\P171\R078\';
         dates = dir(base); dates([1 2])=[];
         Dove_image_file= strcat(base,filesep, dates.name, filesep, '20181204_083439_36_1047_3B_AnalyticMS.tif');
     
     elseif location == 8 % J_SA4
         base='Z:\SpecialNeeds\Sakib\DoveData\P171\R079\';
         dates = dir(base); dates([1 2])=[];
         Dove_image_file= strcat(base,filesep, dates.name, filesep, '20181204_083417_58_1047_3B_AnalyticMS.tif');
     
     elseif location == 9 % CT_SA1
         base='Z:\SpecialNeeds\Sakib\DoveData\P174\R084\';
         dates = dir(base); dates([1 2])=[]; 
         Dove_image_file= strcat(base,filesep, dates.name, filesep, '20180920_092057_1047_3B_AnalyticMS.tif');
     
     elseif location == 10 % CT_SA2
         base='Z:\SpecialNeeds\Sakib\DoveData\P175\R083\';
         dates = dir(base); dates([1 2])=[];
         Dove_image_file= strcat(base,filesep, dates.name, filesep, '20180911_092543_1047_3B_AnalyticMS.tif');

     elseif location == 11 % Auckland_NZ
         base='Z:\SpecialNeeds\Sakib\DoveData\P072\R087\';
         dates = dir(base); dates([1 2])=[];
         Dove_image_file= strcat(base,filesep, dates.name, filesep, '20180909_225853_1047_3B_AnalyticMS.tif');
         
     elseif location == 12 % Sydney_Australia
         base='Z:\SpecialNeeds\Sakib\DoveData\P089\R083\';
         dates = dir(base); dates([1 2])=[];
         Dove_image_file= strcat(base,filesep, dates.name, filesep, '20181120_002253_32_1047_3B_AnalyticMS.tif');
     
     elseif location == 13 % NSW_Aus1
         base='Z:\SpecialNeeds\Sakib\DoveData\P090\R082\';
         dates = dir(base);dates([1 2])=[];
         Dove_image_file= strcat(base,filesep, dates.name, filesep, '20181111_002927_91_1047_3B_AnalyticMS.tif');

         MData_file=strcat(base, filesep, dates.name, filesep, '20181111_002927_91_1047_3B_AnalyticMS_metadata.xml');
         [MData_values]= xml2struct_new_v(MData_file);

     elseif location == 14 % NSW_Aus2
         base='Z:\SpecialNeeds\Sakib\DoveData\P091\R082\';
         dates = dir(base); dates([1 2])=[];
         Dove_image_file= strcat(base,filesep, dates.name, filesep, '20181102_003531_56_1047_3B_AnalyticMS.tif');

     elseif location == 15 % Melbourne_Australia
         base='Z:\SpecialNeeds\Sakib\DoveData\P093\R086\';
         dates = dir(base); dates([1 2])=[];
         Dove_image_file= strcat(base,filesep, dates.name, filesep, '20181116_005548_35_1047_3B_AnalyticMS.tif');
         
     elseif location == 16 % Libya4_S52
         base='Z:\SpecialNeeds\Sakib\DoveData\P181\R040\';
         dates = dir(base); dates([1 2])=[];
         Dove_image_file= strcat(base, filesep, dates.name, filesep, '20181007_082530_52_1047_3B_AnalyticMS.tif');
  
     elseif location == 17 % Libya4_S50
         base='Z:\SpecialNeeds\Sakib\DoveData\P181\R040\';
         dates = dir(base);dates([1 2])=[];
         Dove_image_file= strcat(base,filesep, dates.name, filesep, '20181007_082532_50_1047_3B_AnalyticMS.tif');
         
     elseif location == 18 % Illinois
         base='Z:\SpecialNeeds\Sakib\DoveData\P023\R033\';
         dates = dir(base);dates([1 2])=[];
         Dove_image_file= strcat(base,filesep, dates.name, filesep, '20181121_154822_82_1047_3B_AnalyticMS.tif');
         
     elseif location == 19 % Minnesota
         base='Z:\SpecialNeeds\Sakib\DoveData\P029\R028\';
         dates = dir(base);dates([1 2])=[];
         Dove_image_file= strcat(base,filesep, dates.name, filesep, '20181217_160118_89_1047_3B_AnalyticMS.tif');
         
    elseif location == 20 % South Dakota
         base='Z:\SpecialNeeds\Sakib\DoveData\P029\R029\';
         dates = dir(base);dates([1 2])=[];
         Dove_image_file= strcat(base,filesep, dates.name, filesep, '20181217_160105_03_1047_3B_AnalyticMS.tif');
         
     end
  
  % Reading the RGB-N image
  Dove_image_all_band =imread(Dove_image_file);
  % window = 41; window_L8 = 3;
%%
  for window = 3:2:61
    %for D_band=1:4
        D_band = 1; 
%         % Landsat 8
%         L8_band=D_band+1;
%         Folder_info=dir(fullfile(base, dates.name, strcat('*B', bands{L8_band},'.tif')));  
%         L8_image_file= fullfile(Folder_info.folder, Folder_info.name);
%         
%         % Landsat 8- DN Space
%         L8_image = imread(L8_image_file); L8_info = geotiffinfo(L8_image_file);
%         L8_image = double(L8_image);
%         L8_image(L8_image == 0)= nan;
%         
%         % TOA Radiance
%         MTL=dir(fullfile(base, dates.name,'*MTL.txt'));
%         [MTL_List_L8, vaule]= MTL_parser_L8(fullfile(MTL.folder, MTL.name));
%         
%         rmb= MTL_List_L8.RADIOMETRIC_RESCALING.(strcat('RADIANCE_MULT_BAND_', bands{L8_band}));
%         rab= MTL_List_L8.RADIOMETRIC_RESCALING.(strcat('RADIANCE_ADD_BAND_', bands{L8_band}));
%         TOA_rad_L8 = rmb.*L8_image + rab;
%         mean_TOA = nanmean(nanmean(TOA_rad_L8));
        
        % Dove 
        % Reading the image band-by-band 
        Dove_image=Dove_image_all_band(:,:,D_band);
        Dove_image=double(Dove_image);
        Dove_image(Dove_image == 0)= nan;
        
        band_sf=0.01; % radiometric scale factor
        TOA_rad_D = band_sf.*Dove_image;
        mean_TOA = nanmean(nanmean(TOA_rad_D));
        
        % Reflectance
        %%% Reflectance conversion coefficient
%         if D_band == 1
%             id = 12; 
%         elseif D_band == 2
%             id = 14;
%         elseif D_band == 3
%             id = 16;
%         elseif D_band == 4
%             id = 18;
%         end
%         
%         Ref_coeff =str2num(MData_values_old.Children(10).Children(2).Children(id).Children(10).Children.Data);
%         TOA_ref_D = Dove_image.*Ref_coeff;
        
%         % Dove image as a ROI for Landsat 8
%         UL_x=235772;
%         UL_y=4202871;
%         UR_x=260068;
%         UR_y=4207186;
%         LR_x=262907;
%         LR_y=4190775;
%         LL_x=238671;
%         LL_y=4186460;
% 
%         x_vec=[UL_x UR_x LR_x LL_x UL_x];
%         y_vec=[UL_y UR_y LR_y LL_y UL_y];
% 
%         %Map coordinate to pixel coordinate
%         [Pixel_Row_unrounded, Pixel_Column_unrounded] = map2pix(L8_info.RefMatrix, x_vec, y_vec);
%         Pixel_Row= round(Pixel_Row_unrounded);
%         Pixel_Column= round(Pixel_Column_unrounded);
%         [Row, Column]= size(TOA_rad_L8);
% 
%         % ROI Mask 
%         mask= poly2mask(Pixel_Column, Pixel_Row,  Row, Column);
%         TOA_rad_L8=TOA_rad_L8.*mask;
    
        % Mean and Standard Deviation Image
        Dove_mean_image = conv2(TOA_rad_D, ones(window)./(window*window),'same'); 
        %Dove_mean_image_ref = conv2(TOA_ref_D, ones(window)./(window*window),'same'); 
        %L8_mean_image = conv2(TOA_rad_L8, ones(window_L8)./(window_L8*window_L8),'same'); 
        %mean_TOA = nanmean(nanmean(Dove_mean_image));
        
        % window = 2*k+1; k is defined in this way in movingstd2 
        k= (window-1)/2; 
        %k_L8= (window_L8-1)/2; 
        Dove_SD_image = movingstd2(TOA_rad_D, k);
        %Dove_SD_image_ref = movingstd2(TOA_ref_D, k);
        %L8_SD_image = movingstd2(TOA_rad_L8, k_L8);
        Mean_SD(window) = nanmean(nanmean(Dove_SD_image));
        
        % SNR calculation
        SNR_esti = Dove_mean_image./Dove_SD_image; 
        %SNR_esti_ref = Dove_mean_image_ref./Dove_SD_image_ref; 
        %SNR_esti_L8 = L8_mean_image./L8_SD_image;
        
        % Checking the SNR image and SNR distribution
        %figure; histogram(SNR_esti)
        %figure; histogram(Dove_SD_image)
        
        [N,edges] = histcounts(Dove_SD_image);
        maximum= max(max(N));
        [x, y]= find (N== maximum);
        Noise(window)= edges(y);
       % max_count(window) = max(N);
        window_size(window)= window;
        
%         Noise(Noise==0)= NaN;
%         
%         N(N<=200)= NaN;
%         minimum= min(min(N));
%         [i, j]=find(N==minimum);
%         Noise_highest(window) = edges(max(j));
%         Noise_highest(Noise_highest==0)= NaN;

        
        % subplot(2,2,D_band); % For plotting all bands in a single chart
        % plotting SNR_esti vs. mean TOA Radiance
        %figure, plot(Dove_mean_image, SNR_esti,  '.' ,'color', band_colors{D_band});
%         figure, plot(Dove_mean_image_ref, SNR_esti_ref,  '.' ,'color', band_colors{D_band});
%         figure, plot(L8_mean_image, SNR_esti_L8,  '.' ,'color', band_colors{D_band});
%         
%         % For titles and level
%         title(strcat(strcat(band_name{D_band}, ' Band')));
%         xlabel('Local Standard Deviation')
%         ylabel('Local Mean')
% % %         xlim([0 140]); 
% % %         ylim([0 60])
        hold on
        grid on
        grid minor
        ax  = gca;
        ax.FontSize = 24;
        ax.GridColor = 'k';
%         ax.GridAlpha = 0.8;
        ax.MinorGridColor = 'k';
        ax.FontName = 'Times New Roman';
% %        hold on
   end  
      
    % hold on % for plotting all the bands in a single chart
  end 
%end 
sz= 19;
Noise(Noise==0)= NaN;
figure, plot(window_size, Noise, '.','MarkerSize',sz)
figure, plot(window_size, Noise_highest, '.','MarkerSize',sz)
figure, plot(window_size, Mean_SD, '.','MarkerSize',sz)

leg_location ={'California', 'Wisconsin', 'Indiana',  'Missouri', 'Illinois', 'Minnesota', 'J,SA(P/R-170/078)', 'J,SA(P/R-170/079)', ...
                'J,SA(P/R-171/078)', 'J,SA(P/R-171/079)','CT1,SA(P/R-174/084)', 'CT2,SA(P/R-175/083)',...
               'Auckland,NZ', 'Sydney, Australia', 'NSW,Aus(P/R-090/082)', 'NSW,Aus(P/R-091/082)',...
               'Melbourne, Australia', 'Libya4'};
marker = {'o', 's', 'd', '+','p',  '^','>', '<', 'v', 's', 'd', '+','p','x','*','h','.', 'o','+'};
size= 49;
SNR_all_scene = [81, 62.86, 62.05, 41, 51.25, 55.36, 81.25, 44.27, 75.9, 56.07, 74.53,  72.75, 68, 41.08, 35.28, 64.76, 52.7, 121];
TOArad_all_scene = [47, 44, 45.92,  39.42, 41, 31,  74.75,  87.65, 75.9, 86.35, 52.17, 52.38, 45, 76.88, 74.08, 77.71, 60.8, 123];
figure
for leg =1:18
    scatter(TOArad_all_scene(leg), SNR_all_scene(leg),size, marker{leg},'LineWidth',1)
    hold on
end 

legend(leg_location,'FontSize',16, 'Location','southeast');
ylim([0 150])
xlim([0 130])

hold on
title('Average SNR vs. Mean TOA Radiance')
ylabel ('Average SNR')
xlabel ('Mean TOA Radiance (W/Sr/m^2/{\mum})')

[p_slope, S] = polyfit(Scene, SNR_all_scene, 1);
p_slope = round(p_slope, 5);
fitline = fit(Scene, SNR_all_scene, 'poly1');
f = polyval(p_slope, Scene); plot(Scene, f ); hold on
hold on        
title('Average SNR of 1047 Dove Sensor')
ylabel('Average SNR')
xlabel('Different Scenes')

hold on; ylim([20 120]); ylim([0 1.2])
title('Local SNR vs. Local Mean TOA Reflectance')
ylabel('Local SNR')
%xlabel('Local Mean TOA Radiance (W/Sr/m^2/{\mum})')
xlabel('Local Mean TOA Reflectance')

hold on
xlim([0 6])
title('Noise Distribution')
xlabel('Standard Deviation')
ylabel('No. of Observation')


% figure, imshow(SNR_esti, [])
% figure, imagesc(SNR_esti, [0 1000]); colorbar

hold on
title('Mean Noise vs. Window Size')
ylabel ('Mean Noise')
xlabel ('Window Size')
ylim([0 1.3]);


filename = 'Z:\SpecialNeeds\Sakib\DoveData\April Data\2019_Q1_R2_QUALITY_DATA_PACKET_scenes.json';
val = jsondecode(fileread(filename));
