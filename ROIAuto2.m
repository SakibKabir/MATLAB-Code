clear all
%%% Directory 
for location = 1:19
     if location == 1 % California
         base='Z:\ImageDrive\PlanetLabs\Processed\1047\P039\R037\';  
     elseif location == 2 % Wisconsin
         base='Z:\ImageDrive\PlanetLabs\Processed\1047\P026\R029\';  
     elseif location == 3 % Indiana
         base='Z:\ImageDrive\PlanetLabs\Processed\1047\P021\R032\';
     elseif location == 4 % Missouri
         base='Z:\ImageDrive\PlanetLabs\Processed\1047\P023\R034\'; 
     elseif location == 5 % Illinois
         base='Z:\ImageDrive\PlanetLabs\Processed\1047\P023\R033\';
     elseif location == 6 % Minnesota
         base='Z:\ImageDrive\PlanetLabs\Processed\1047\P029\R028\';  
     elseif location == 7 % South Dakota
         base='Z:\ImageDrive\PlanetLabs\Processed\1047\P029\R029\';      
     elseif location == 8 % J_SA1
         base='Z:\ImageDrive\PlanetLabs\Processed\1047\P170\R078\';
     elseif location == 9 % J_SA2
         base='Z:\ImageDrive\PlanetLabs\Processed\1047\P170\R079\';
     elseif location == 10 % J_SA3
         base='Z:\ImageDrive\PlanetLabs\Processed\1047\P171\R078\';
     elseif location == 11 % J_SA4
         base='Z:\ImageDrive\PlanetLabs\Processed\1047\P171\R079\';
     elseif location == 12 % CT_SA1
         base='Z:\ImageDrive\PlanetLabs\Processed\1047\P174\R084\';     
     elseif location == 13 % CT_SA2
         base='Z:\ImageDrive\PlanetLabs\Processed\1047\P175\R083\';
     elseif location == 14 % Auckland_NZ
         base='Z:\ImageDrive\PlanetLabs\Processed\1047\P072\R087\';
     elseif location == 15 % Sydney_Australia
         base='Z:\ImageDrive\PlanetLabs\Processed\1047\P089\R083\';
     elseif location == 16 % NSW_Aus1
         base='Z:\ImageDrive\PlanetLabs\Processed\1047\P090\R082\';
     elseif location == 17 % NSW_Aus2
         base='Z:\ImageDrive\PlanetLabs\Processed\1047\P091\R082\';
     elseif location == 18 % Melbourne_Australia
         base='Z:\ImageDrive\PlanetLabs\Processed\1047\P093\R086\';
     elseif location == 19 % Libya4
         base='Z:\ImageDrive\PlanetLabs\Processed\1047\P181\R040\';
     end
%%
% clear all
dates = dir(base); dates([1 2])=[];

  %%% SNR at Typical Radiance level at +/- 1.25%
  Ltyp = 1;
  %%% SNR at Typical Radiance level at +/- 2.5%   
  % Ltyp = 2;

% Dove File Path
DoveDataPath = fullfile(base, dates.name, '*_AnalyticMS.tif');
DoveDirectory = dir(DoveDataPath);
% DoveDataPathXML = fullfile(base, dates.name, '*_metadata.xml');
DoveDataPathXML = fullfile(base, dates.name, '*_metadata.xml');
DoveDirectoryXML = dir(DoveDataPathXML);
[NoOf_Scene, y] = size(DoveDirectory);

%% Dove Scenes
for D_band = 1 %:4
 for scene = 1 %: NoOf_Scene
    % Extracting all the values from Metadata file
    MData_file = fullfile(base, dates.name, DoveDirectoryXML(scene).name);
    [MData_values]= xml2struct_new_v(MData_file);
    [MData_values_old]= xml2struct(MData_file);
    
    %%% Dove
    % Dove Image file- scene by scene
    Dove_image_file = fullfile(base, dates.name, DoveDirectory(scene).name);

    % Reading Dove Image
    Dove_image_all_band =imread(Dove_image_file);
    DN_Dove = Dove_image_all_band(:,:, D_band);
    DN_Dove = double(DN_Dove);
    DN_Dove(DN_Dove == 0) = NaN;
    band_sf = 0.01;
    Dove_TOArad = DN_Dove*band_sf;
    
    % figure, imagesc(Dove_TOArad); colorbar
    % figure, histogram(Dove_TOArad)

    %% Dove Image Processing
    %%%================== Local Mean Section
    window = 29;
    k= (window-1)/2; 
    Dove_mean_image = conv2(Dove_TOArad, ones(window)./(window*window),'same');
    % figure, imagesc(Dove_mean_image); colorbar
    % figure, histogram(Dove_mean_image)
    Dove_mean_image_filtered = Dove_mean_image;
    
    %%% SNR at Typical Radiance level at +/- 1.25%
    if Ltyp == 1
        if D_band == 1
            L_low = 39.5; L_high = 40.5; % Blue - 40:
        elseif D_band == 2 
            L_low = 29.6; L_high = 30.4; % Green - 30:
        elseif D_band == 3
            L_low = 21.7; L_high = 22.3; % Red - 22:
        elseif D_band == 4
            L_low = 13.8; L_high = 14.2; % NIR - 14:
        end 
        
    %%% SNR at Typical Radiance level at +/- 2.5%
    elseif Ltyp == 2 
        if D_band == 1
            L_low = 39; L_high = 41; % Blue - 40:
        elseif D_band == 2 
            L_low = 29.3; L_high = 30.8; % Green - 30:
        elseif D_band == 3
            L_low = 21.5; L_high = 22.6; % Red - 22:
        elseif D_band == 4
            L_low = 13.7; L_high = 14.4; % NIR - 14:
        end 
    end
    
    % Dove_mean_image_filtered(Dove_mean_image_filtered < L_low | Dove_mean_image_filtered > L_high) = NaN;
    Dove_TOArad_fil = Dove_TOArad;
    Dove_TOArad_fil(Dove_TOArad_fil < L_low | Dove_TOArad_fil > L_high) = NaN;
    % figure, imagesc(Dove_TOArad_fil); colorbar
    
    %%%=================== Local SD Section
    D_LSD_image = movingstd2(Dove_TOArad, k);
    % figure, imagesc(D_LSD_image); colorbar
    % figure, histogram(D_LSD_image);
    
    %%%================ SNR Section
    % D_SNR_filtered = Dove_mean_image_filtered./D_LSD_image;
     D_SNR_filtered = Dove_TOArad_fil./D_LSD_image;
    % figure, imagesc(D_SNR_filtered); colorbar
    % figure, histogram(D_SNR_filtered)
    Data_prc = sum(sum(~isnan(D_SNR_filtered)))/sum(sum(~isnan(Dove_mean_image)))*100;
    
    %%% Mean, percentage calculation
    temp_D_SNR_filtered = D_SNR_filtered;
    temp_D_SNR_filtered = temp_D_SNR_filtered(~isnan(temp_D_SNR_filtered));
    Mean = round(mean(temp_D_SNR_filtered), 1);
    Std_D = std(temp_D_SNR_filtered);
    Median = median(temp_D_SNR_filtered);
    
    % Percentage calculation
    temp_D_SNR_filtered_sorted = sort(temp_D_SNR_filtered);
    ind_95 = round(0.95*length(temp_D_SNR_filtered_sorted));
    SNR_95prc = temp_D_SNR_filtered_sorted(find(temp_D_SNR_filtered_sorted) == ind_95);
    
    ind_96 = round(0.96*length(temp_D_SNR_filtered_sorted));
    SNR_96prc = temp_D_SNR_filtered_sorted(find(temp_D_SNR_filtered_sorted) == ind_96);
    
    ind_97 = round(0.97*length(temp_D_SNR_filtered_sorted));
    SNR_97prc = temp_D_SNR_filtered_sorted(find(temp_D_SNR_filtered_sorted) == ind_97);
    
    ind_98 = round(0.98*length(temp_D_SNR_filtered_sorted));
    SNR_98prc = temp_D_SNR_filtered_sorted(find(temp_D_SNR_filtered_sorted) == ind_98);
    
    ind_99 = round(0.99*length(temp_D_SNR_filtered_sorted));
    SNR_99prc = temp_D_SNR_filtered_sorted(find(temp_D_SNR_filtered_sorted) == ind_99);
    
    ind_99_5 = round(0.995*length(temp_D_SNR_filtered_sorted));
    SNR_99_5prc = temp_D_SNR_filtered_sorted(find(temp_D_SNR_filtered_sorted) == ind_99_5);
    
    ind_99_7 = round(0.997*length(temp_D_SNR_filtered_sorted));
    SNR_99_7prc = temp_D_SNR_filtered_sorted(find(temp_D_SNR_filtered_sorted) == ind_99_7);
    
    ind_99_8 = round(0.998*length(temp_D_SNR_filtered_sorted));
    SNR_99_8prc = temp_D_SNR_filtered_sorted(find(temp_D_SNR_filtered_sorted) == ind_99_8);
    
    ind_99_9 = round(0.999*length(temp_D_SNR_filtered_sorted));
    SNR_99_9prc = temp_D_SNR_filtered_sorted(find(temp_D_SNR_filtered_sorted) == ind_99_9);
    
%     SDs = round((Mean + (1:3)*Std_D),1);

     %% Plotting section
%     figure(scene)
%     h = histogram(D_SNR_filtered); % Histogram
%     N = max(h.Values); % Maximum bin count
% 
%     hold on
%     plot([Mean Mean],[0 N],'r','LineWidth',2) % 97 percentage
%     plot([SNR_98prc SNR_98prc],[0 N/3],'g','LineWidth',3) % 98 percentage
%     plot([SNR_99prc SNR_99prc],[0 N/3],'m','LineWidth',3) % 99 percentage
%     plot([SNR_99_5prc SNR_99_5prc],[0 N/3],'c','LineWidth',3) % 99.5 percentage
% %     plot([SNR_99_7prc SNR_99_7prc],[0 N/4],'b','LineWidth',3) % 99.7 percentage
% %     plot([SNR_99_8prc SNR_99_8prc],[0 N/4],'g','LineWidth',3) % 99.8 percentage
% %     plot([SNR_99_9prc SNR_99_9prc],[0 N/4],'k','LineWidth',3) % 99.9 percentage
%     xlim([0 120])
%     hold on
%     % Legend for L_typical
%     legend([num2str(Data_prc), '% of total data'], ['Mean (' num2str(Mean) ')'],...
%         ['SNR at 98% (' num2str(SNR_98prc) ')'],...
%         ['SNR at 99% (' num2str(SNR_99prc) ')'], ['SNR at 99.5% (' num2str(SNR_99_5prc) ')']);
%     
%     % Legend for L_high
% %     legend([num2str(Data_prc), '% of total data'], ['Mean (' num2str(Mean) ')'],...
% %         ['SNR at 99% (' num2str(SNR_99prc) ')'],...
% %         ['SNR at 99.5% (' num2str(SNR_99_5prc) ')'], ['SNR at 99.7% (' num2str(SNR_99_7prc) ')'],...
% %         ['SNR at 99.8% (' num2str(SNR_99_8prc) ')'], ['SNR at 99.9% (' num2str(SNR_99_9prc) ')']);
%     
%     title('SNR Distribution')
%     xlabel('Signal to Noise Ratio')
%     ylabel('Number of Observation')
%     ax  = gca;
%     grid on
%     ax.FontSize = 36;
%     hold off
    
     %% Storing the data
     if Data_prc > 0.5
         SNR_95(scene) = SNR_95prc;
         SNR_96(scene) = SNR_96prc;
         SNR_97(scene) = SNR_97prc;
         SNR_98(scene) = SNR_98prc;
         SNR_99(scene) = SNR_99prc;
         SNR_99_5(scene) = SNR_99_5prc;
         SNR_99_7(scene) = SNR_99_7prc;
         SNR_99_8(scene) = SNR_99_8prc;
         SNR_99_9(scene) = SNR_99_9prc;
     else
         SNR_95(scene) = NaN;
         SNR_96(scene) = NaN;
         SNR_97(scene) = NaN;
         SNR_98(scene) = NaN;
         SNR_99(scene) = NaN;
         SNR_99_5(scene) = NaN;
         SNR_99_7(scene) = NaN;
         SNR_99_8(scene) = NaN;
         SNR_99_9(scene) = NaN;
     end

%     Per_Data(scene) = Data_prc;
%     SNR_all_99prc(scene) = SNR_99prc;

end
     DoveSNR.SNR_95_D1047{location,D_band} = SNR_95;
     DoveSNR.SNR_96_D1047{location,D_band} = SNR_96;
     DoveSNR.SNR_97_D1047{location,D_band} = SNR_97;
     DoveSNR.SNR_98_D1047{location,D_band} = SNR_98;
     DoveSNR.SNR_99_D1047{location,D_band} = SNR_99;
     DoveSNR.SNR_99_5_D1047{location,D_band} = SNR_99_5;
     DoveSNR.SNR_99_7_D1047{location,D_band} = SNR_99_7;
     DoveSNR.SNR_99_8_D1047{location,D_band} = SNR_99_8;
     DoveSNR.SNR_99_9_D1047{location,D_band} = SNR_99_9;
end
clearvars -except DoveSNR
end
%% 
% load('DoveSNR_B_1prc.mat')
SNRall_97_mat = nan(19, 10);
SNRall_98_mat = nan(19, 10);
SNRall_99_mat = nan(19, 10);
SNRall_99_5_mat = nan(19, 10);
SNRall_99_7_mat = nan(19, 10);
SNRall_99_8_mat = nan(19, 10);
SNRall_99_9_mat = nan(19, 10);

for location = 1:19
    SNRall_97 = DoveSNR.SNR_97_D1047{location, 1};
    SNRall_97_mat(location, 1:length(SNRall_97)) = SNRall_97;
    
    SNRall_98 = DoveSNR.SNR_98_D1047{location, 1};
    SNRall_98_mat(location, 1:length(SNRall_98)) = SNRall_98;
    
    SNRall_99 = DoveSNR.SNR_99_D1047{location, 1};
    SNRall_99_mat(location, 1:length(SNRall_99)) = SNRall_99;
    
    SNRall_99_5 = DoveSNR.SNR_99_5_D1047{location, 1};
    SNRall_99_5_mat(location, 1:length(SNRall_99_5)) = SNRall_99_5;
    
    SNRall_99_7 = DoveSNR.SNR_99_7_D1047{location, 1};
    SNRall_99_7_mat(location, 1:length(SNRall_99_7)) = SNRall_99_7;
    
    SNRall_99_8 = DoveSNR.SNR_99_8_D1047{location, 1};
    SNRall_99_8_mat(location, 1:length(SNRall_99_8)) = SNRall_99_8;
        
    SNRall_99_9 = DoveSNR.SNR_99_9_D1047{location, 1};
    SNRall_99_9_mat(location, 1:length(SNRall_99_9)) = SNRall_99_9;
end

%% Box Plot
figure,
b = boxplot([SNRall_97_mat(:), SNRall_98_mat(:), SNRall_99_mat(:), SNRall_99_5_mat(:)],...
         'Labels', {'97%','98%','99%', '99.5%'},'Symbol', 'b.', 'OutlierSize',30,'Whisker',1.5);
set(b,{'linew'},{2})
%ylim([0 100])
xlabel('SNR at Percentage')
ylabel('Signal to Noise Ratio')
title('Signal to Noise Ratio of Dove 1047 NIR Band')

ax = gca;
ax.FontSize = 36;
ylim([0 100])
grid on;
grid minor
ax.GridColor = 'k';
ax.MinorGridColor = 'k';
     
%% Plot
    figure,
    range = ones(1, length(SNRall_98_mat(:)));
    x = 0.4*range;
    plot(x, SNRall_98_mat(:), 'b.', 'MarkerSize', 30, 'LineWidth', 3);
    hold on;
    x = 0.6*range;
    plot(x, SNRall_99_mat(:), 'g.', 'MarkerSize', 30, 'LineWidth', 3);
    hold on; 
    x = 0.8*range;
    plot(x, SNRall_99_5_mat(:), 'm.', 'MarkerSize', 30, 'LineWidth', 3);
    
    % Set up axes.
    xlim([0, 1.2]);
    ylim([0, 120]);
    xlabel('SNR at Percentage')
    ylabel('Signal to Noise Ratio')
    title('Signal to Noise Ratio of Dove 1047')
    
    ax = gca;
    ax.XTick = [0.4, 0.6, 0.8];
    ax.XTickLabels = {'98%', '99%', '99.5%'};
    ax.FontSize = 36;
    grid on;
    grid minor
    ax.GridColor = 'k';
%   ax.GridAlpha = 0.8;
    ax.MinorGridColor = 'k';
    %ax.FontName = 'Times New Roman';
    
%% Dove SNR all band together
for band = 1:3
       if band == 1
            load('DoveSNR_B_1prc.mat')
            % Creating Nan matrix to store the values
            SNRall_97_mat = nan(19, 10);
            SNRall_98_mat = nan(19, 10); SNRall_99_mat = nan(19, 10);
            SNRall_99_5_mat = nan(19, 10); SNRall_99_7_mat = nan(19, 10);
        elseif band == 2
            load('DoveSNR_G_1prc.mat')
            % Creating Nan matrix to store the values
            SNRall_97_mat = nan(19, 10);
            SNRall_98_mat = nan(19, 10); SNRall_99_mat = nan(19, 10);
            SNRall_99_5_mat = nan(19, 10); SNRall_99_7_mat = nan(19, 10);      
        elseif band == 3
            load('DoveSNR_R_1prc.mat')
            % Creating Nan matrix to store the values
            SNRall_97_mat = nan(19, 10);
            SNRall_98_mat = nan(19, 10); SNRall_99_mat = nan(19, 10);
            SNRall_99_5_mat = nan(19, 10); SNRall_99_7_mat = nan(19, 10);
        elseif band == 4
            load('DoveSNR_NIR.mat')
            % Creating Nan matrix to store the values
            SNRall_97_mat = nan(19, 10);
            SNRall_98_mat = nan(19, 10); SNRall_99_mat = nan(19, 10);
            SNRall_99_5_mat = nan(19, 10); SNRall_99_7_mat = nan(19, 10);  
        end
    
    for location = 1:19
        SNRall_97 = DoveSNR.SNR_97_D1047{location, 1};
        SNRall_97_mat(location, 1:length(SNRall_97)) = SNRall_97;
        
        SNRall_98 = DoveSNR.SNR_98_D1047{location, 1};
        SNRall_98_mat(location, 1:length(SNRall_98)) = SNRall_98;

        SNRall_99 = DoveSNR.SNR_99_D1047{location, 1};
        SNRall_99_mat(location, 1:length(SNRall_99)) = SNRall_99;

        SNRall_99_5 = DoveSNR.SNR_99_5_D1047{location, 1};
        SNRall_99_5_mat(location, 1:length(SNRall_99_5)) = SNRall_99_5;
        
        SNRall_99_7 = DoveSNR.SNR_99_7_D1047{location, 1};
        SNRall_99_7_mat(location, 1:length(SNRall_99_7)) = SNRall_99_7;
    end
    
        SNRall_97_vec(:, band) = SNRall_97_mat(:);
        SNRall_98_vec(:, band) = SNRall_98_mat(:);
        SNRall_99_vec(:, band) = SNRall_99_mat(:);
        SNRall_99_5_vec(:, band) = SNRall_99_5_mat(:);
        SNRall_99_7_vec(:, band) = SNRall_99_7_mat(:);
 end

%% Ploting all band together
band_name = {'BLUE','GREEN' ,'RED' ,'NIR' };

for band = 1:3
    figure(band)
    b = boxplot([SNRall_98_vec(:,band), SNRall_99_vec(:,band), SNRall_99_5_vec(:,band), SNRall_99_7_vec(:,band)],...
         'Labels', {'98%','99%', '99.5%', '99.7%' },'Symbol','b.',  'OutlierSize',30,'Whisker',1.5);
    hold on
    ylim([0 100])
    set(b,{'linew'},{2})
    xlabel('SNR at Percentage')
    ylabel('Signal to Noise Ratio')
    title(strcat('Signal to Noise Ratio of Dove 1047',  {', '}, strcat(band_name{band}, ' Band')))
    ax  = gca;
    grid on
    grid minor
    ax.FontSize = 36;
    hold off
end

%% Ploting Single SNR
    figure(1)
    b = boxplot([SNRall_97_vec(:,1), SNRall_99_5_vec(:,2), SNRall_99_8_vec(:,3)],...
         'Labels', {'Blue', 'Green', 'Red'},'Symbol','b.',  'OutlierSize',30,'Whisker',1.5);
    blue_q3 = quantile(SNRall_98_vec(:,1), 0.75);
    green_q3 = quantile(SNRall_98_vec(:,1), 0.75);
    hold on
    ylim([0 100])
    set(b,{'linew'},{2})
    %xlabel('BanSNR')
    ylabel('Signal to Noise Ratio')
    title('Signal to Noise Ratio of Dove 1047')
    ax  = gca;
    grid on
    grid minor
    ax.FontSize = 36;
    hold off
