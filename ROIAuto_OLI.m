clear all
for location = 1:17
    if location == 1
        base='Z:\ImageDrive\OLI.TIRS\L8\P023\R032';
    elseif location == 2
        base='Z:\ImageDrive\OLI.TIRS\L8\P023\R033';
    elseif location == 3
        base='Z:\ImageDrive\OLI.TIRS\L8\P023\R034';
    elseif location == 4
        base='Z:\ImageDrive\OLI.TIRS\L8\P023\R035';
    elseif location == 5
        base='Z:\ImageDrive\OLI.TIRS\L8\P023\R036';
    elseif location == 6
        base='Z:\ImageDrive\OLI.TIRS\L8\P023\R037';
    elseif location == 7    
        base='Z:\ImageDrive\OLI.TIRS\L8\P026\R029';
    elseif location == 8    
        base='Z:\ImageDrive\OLI.TIRS\L8\P021\R032';
    elseif location == 9    
         base='Z:\ImageDrive\OLI.TIRS\L8\P029\R028';
    elseif location == 10    
         base='Z:\ImageDrive\OLI.TIRS\L8\P022\R039';
    elseif location == 11    
         base='Z:\ImageDrive\OLI.TIRS\L8\P072\R087';
    elseif location == 12    
         base='Z:\ImageDrive\OLI.TIRS\L8\P089\R083';
    elseif location == 13    
         base='Z:\ImageDrive\OLI.TIRS\L8\P090\R082';     
    elseif location == 14    
         base='Z:\ImageDrive\OLI.TIRS\L8\P091\R082';    
    elseif location == 15    
         base='Z:\ImageDrive\OLI.TIRS\L8\P093\R086';
    elseif location == 16    
         base='Z:\ImageDrive\OLI.TIRS\L8\P171\R078';
    elseif location == 17    
         base='Z:\ImageDrive\OLI.TIRS\L8\P171\R079';   
    elseif location == 18    
         base='Z:\ImageDrive\OLI.TIRS\L8\P038\R037';
    elseif location == 19    
         base='Z:\ImageDrive\OLI.TIRS\L8\P029\R029';
    elseif location == 20    
         base='Z:\ImageDrive\OLI.TIRS\L8\P039\R037';
    elseif location == 21    
         base='Z:\ImageDrive\OLI.TIRS\L8\P001\R065';
    elseif location == 22
         base='Z:\ImageDrive\OLI.TIRS\L8\P023\R033';
    elseif location == 23 % Libya 4
         base='Z:\ImageDrive\OLI.TIRS\L8\P181\R040';
    elseif location == 24 % Niger 1
         base='Z:\ImageDrive\OLI.TIRS\L8\P189\R046';    
    elseif location == 25 % Niger 2
         base='Z:\ImageDrive\OLI.TIRS\L8\P188\R045'; 
    elseif location == 26 % Sonoran Dessert
         base='Z:\ImageDrive\OLI.TIRS\L8\P038\R038'; 
    elseif location == 27 %
         base='Z:\ImageDrive\OLI.TIRS\L8\P190\R043';
    elseif location % new locations    
         base='Z:\ImageDrive\OLI.TIRS\L8\P022\R032';
    elseif location % new locations    
         base='Z:\ImageDrive\OLI.TIRS\L8\P022\R033';
    end
    
  %%% clear all
  %%% SNR at Typical Radiance level at +/- 1.25%
   Ltyp = 1;
  %%% SNR at Typical Radiance level at +/- 2.5%   
  % Ltyp = 2;
  dates = dir(base); dates([1 2])=[]; 
  bands =  {'1' '2' '3' '4' '5' }; 
   
for band = 5:5
 for date = 1: size(dates,1)
  if exist (fullfile(base, dates(date).name, 'LC1'),'dir')
    %% MTL File Parsing 
    MTL = dir(fullfile(base, dates(date).name,'LC1','*MTL.txt'));
    [MTL_List, vaule]= MTL_parser_L8(fullfile(MTL.folder, MTL.name));

    %Image file location
    Band_No = dir(fullfile(base, dates(date).name, 'LC1', strcat('*B', bands{band},'.tif')));
    L8_image_file = strcat(base, filesep, dates(date).name, filesep, 'LC1', filesep, Band_No.name);
    
    % BQA File location and Reading the BQA Image, creating mask
    BQA_File_dir = dir(fullfile(base, dates(date).name, 'LC1', '*_BQA.TIF'));
    L8_BQA_File = strcat(base, filesep, dates(date).name, filesep, 'LC1', filesep, BQA_File_dir.name);
    [BQA_data, R] = geotiffread(L8_BQA_File);
    BQA_mask =(BQA_data == 2720);
    
    % Reading the OLI image
    [DN_L8, R_L8] = geotiffread(L8_image_file);
    DN_L8 = double(DN_L8);
    DN_L8(DN_L8 == 0)= nan;
    L8_image_info = geotiffinfo(L8_image_file);
    % figure, imagesc(DN_L8); colorbar
    
    % Mulitplicative and additive factors
    rmb = MTL_List.RADIOMETRIC_RESCALING.(strcat('RADIANCE_MULT_BAND_', bands{band}));
    rab = MTL_List.RADIOMETRIC_RESCALING.(strcat('RADIANCE_ADD_BAND_', bands{band}));
    
    % L8 TOA Radiance
    L8_TOArad_raw = DN_L8*rmb+rab;
    % figure, histogram(L8_TOArad_raw); 
    
    % L8 TOA Radiance masked
    L8_TOArad = L8_TOArad_raw.*BQA_mask;
    L8_TOArad(L8_TOArad == 0)= nan;
    % figure, imagesc(L8_TOArad); colorbar
    % figure, histogram(L8_TOArad);
    
   end

%%  %%%% L8 SNR Calculation 
    window = 3;
    k= (window-1)/2; 
    % L8_mean_image = conv2(L8_TOArad, ones(window)./(window*window),'same');
    % figure, imagesc(L8_mean_image); colorbar
    % figure, histogram(L8_mean_image)
    
    %%% Filtering at typical/high radiance level
    %L8_mean_image_filtered = L8_mean_image;
    
    %%% SNR at Typical Radiance level at +/- 1.25%
    if Ltyp == 1
        if band == 2
            L_low = 39.5; L_high = 40.5; % Blue - 40:
        elseif band == 3 
            L_low = 29.6; L_high = 30.4; % Green - 30:
        elseif band == 4
            L_low = 21.7; L_high = 22.3; % Red - 22:
        elseif band == 5
            L_low = 13.8; L_high = 14.2; % NIR - 14:
        end 
    %%% SNR at Typical Radiance level at +/- 2.5%
    elseif Ltyp == 2 
        if band == 2
            L_low = 39; L_high = 41; % Blue - 40:
        elseif band == 3 
            L_low = 29.3; L_high = 30.8; % Green - 30:
        elseif band == 4
            L_low = 21.5; L_high = 22.6; % Red - 22:
        elseif band == 5
            L_low = 13.7; L_high = 14.4; % NIR - 14:
        end 
    end
    
%     %%% SNR at High Radiance level at +/- 1.25%
%     if L == 1
%         if band == 2
%             L_low = 187.6; L_high = 192.4; % Blue - 190:
%         elseif band == 3 
%             L_low = 191.6; L_high = 196.4; % Green - 194:
%         elseif band == 4
%             L_low = 148.1; L_high = 151.9; % Red - 150:
%         elseif band == 5
%             L_low = 148.1; L_high = 151.9; % NIR - 150:
%         end 
%     %%% SNR at High Radiance level at +/- 2.5%
%     else 
%         if band == 2
%             L_low = 185.3; L_high = 194.8; % Blue - 190:
%         elseif band == 3 
%             L_low = 189.2; L_high = 198.9; % Green - 194:
%         elseif band == 4
%             L_low = 146.3; L_high = 153.8; % Red - 150:
%         elseif band == 5
%             L_low = 146.3; L_high = 153.8; % NIR - 150:
%         end 
%     end

    %L8_mean_image_filtered(L8_mean_image_filtered < L_low | L8_mean_image_filtered > L_high) = NaN;
    
     L8_TOArad_fil = L8_TOArad;
     L8_TOArad_fil(L8_TOArad_fil < L_low | L8_TOArad_fil > L_high) = NaN;
    % figure, imagesc(L8_mean_image_filtered); colorbar
    % figure, histogram(L8_mean_image_filtered)

    %%% L8 Standard Deviation image
    L8_SD_image = movingstd2(L8_TOArad, k);
    % figure, imagesc(L8_SD_image); colorbar
    % figure, histogram(L8_SD_image)

    %%% SNR 
    %L8_SNR_filtered = L8_mean_image_filtered./L8_SD_image;
    L8_SNR_filtered = L8_TOArad_fil./L8_SD_image;
    % figure, imagesc(L8_SNR_filtered); colorbar
    
   % figure('units','normalized','outerposition',[0 0 1 1],'visible','off')
   % figure,
   % figure,imagesc(L8_SNR_filtered, [0 400]); colorbar
   % FileNString = strcat('Z:\SpecialNeeds\Sakib\Figures\', strcat(base(end-8:end-5),...
       % base(end-3:end),'d', dates(date).name, 'band', num2str(band)));
   % print(FileNString, '-djpeg','-r500')
   % figure, histogram(L8_SNR_filtered)
    
    Data_prc = sum(sum(~isnan(L8_SNR_filtered)))/sum(sum(~isnan(L8_TOArad)))*100;

    %%% Meand and SD calculation
    temp_L8_SNR_filtered = L8_SNR_filtered;
    temp_L8_SNR_filtered = temp_L8_SNR_filtered(~isnan(temp_L8_SNR_filtered));
    Mean = round(mean(temp_L8_SNR_filtered), 1);

    temp_L8_SNR_filtered_sorted = sort(temp_L8_SNR_filtered);
    ind_95 = round(0.95*length(temp_L8_SNR_filtered_sorted));
    SNR_95prc = temp_L8_SNR_filtered_sorted(find(temp_L8_SNR_filtered_sorted) == ind_95);
    
    ind_96 = round(0.96*length(temp_L8_SNR_filtered_sorted));
    SNR_96prc = temp_L8_SNR_filtered_sorted(find(temp_L8_SNR_filtered_sorted) == ind_96);
    
    temp_L8_SNR_filtered_sorted = sort(temp_L8_SNR_filtered);
    ind_97 = round(0.97*length(temp_L8_SNR_filtered_sorted));
    SNR_97prc = temp_L8_SNR_filtered_sorted(find(temp_L8_SNR_filtered_sorted) == ind_97);
    
    temp_L8_SNR_filtered_sorted = sort(temp_L8_SNR_filtered);
    ind_98 = round(0.98*length(temp_L8_SNR_filtered_sorted));
    SNR_98prc = temp_L8_SNR_filtered_sorted(find(temp_L8_SNR_filtered_sorted) == ind_98);

    ind_99 = round(0.99*length(temp_L8_SNR_filtered_sorted));
    SNR_99prc = temp_L8_SNR_filtered_sorted(find(temp_L8_SNR_filtered_sorted) == ind_99);
    
    ind_99_5 = round(0.995*length(temp_L8_SNR_filtered_sorted));
    SNR_99_5prc = temp_L8_SNR_filtered_sorted(find(temp_L8_SNR_filtered_sorted) == ind_99_5);
   
    ind_99_7 = round(0.997*length(temp_L8_SNR_filtered_sorted));
    SNR_99_7prc = temp_L8_SNR_filtered_sorted(find(temp_L8_SNR_filtered_sorted) == ind_99_7);
    
    ind_99_8 = round(0.998*length(temp_L8_SNR_filtered_sorted));
    SNR_99_8prc = temp_L8_SNR_filtered_sorted(find(temp_L8_SNR_filtered_sorted) == ind_99_8);
    
    ind_99_9 = round(0.999*length(temp_L8_SNR_filtered_sorted));
    SNR_99_9prc = temp_L8_SNR_filtered_sorted(find(temp_L8_SNR_filtered_sorted) == ind_99_9);
   
%%     %%%% Plotting SNR Distribution
%     if nansum(nansum(L8_SNR_filtered)) > 1
%     %figure('units','normalized','outerposition',[0 0 1 1],'visible','off') 
    figure, 
    L8_SNR_pre = {'0', '370', '305', '230', '205'}; % At Typical Radiance Level
    
    h = histogram(L8_SNR_filtered); % Histogram
    N = max(h.Values); % Maximum bin count
    hold on
    plot([Mean Mean],[0 N],'y','LineWidth',2) % 97 percentage
    plot([str2double(L8_SNR_pre{band}) str2double(L8_SNR_pre{band})], [0 N/2], 'r','LineWidth',2.5)
    %plot([SNR_96prc SNR_98prc],[0 N/5],'g','LineWidth',3) % 96 percentage
    plot([SNR_97prc SNR_97prc],[0 N/5],'b','LineWidth',3) % 97 percentage
    plot([SNR_98prc SNR_98prc],[0 N/5],'g','LineWidth',3) % 98 percentage
    plot([SNR_99prc SNR_99prc],[0 N/5],'m','LineWidth',3) % 99 percentage
    plot([SNR_99_5prc SNR_99_5prc],[0 N/5],'c','LineWidth',3) % 99.5 percentage
%     plot([SNR_99_7prc SNR_99_7prc],[0 N/4],'b','LineWidth',3) % 99.7 percentage
%     plot([SNR_99_8prc SNR_99_8prc],[0 N/4],'g','LineWidth',3) % 99.8 percentage
%     plot([SNR_99_9prc SNR_99_9prc],[0 N/4],'k','LineWidth',3) % 99.9 percentage
    hold on
    
    % Legend for L_typical
    legend([num2str(Data_prc), '% of total data'], ['Mean (' num2str(Mean) ')'],...
           ['Pre-Launch SNR (' L8_SNR_pre{band} ')'],...
           ['SNR at 97% (' num2str(SNR_97prc) ')'],['SNR at 98% (' num2str(SNR_98prc) ')'],...
           ['SNR at 99% (' num2str(SNR_99prc) ')'], ['SNR at 99.5% (' num2str(SNR_99_5prc) ')']);

    % Legend for L_high
%   legend([num2str(Data_prc), '% of total data'], ['Mean (' num2str(Mean) ')'],...
%         ['SNR at 99% (' num2str(SNR_99prc) ')'],...
%         ['SNR at 99.5% (' num2str(SNR_99_5prc) ')'], ['SNR at 99.7% (' num2str(SNR_99_7prc) ')'],...
%         ['SNR at 99.8% (' num2str(SNR_99_8prc) ')'], ['SNR at 99.9% (' num2str(SNR_99_9prc) ')']);
    
    % Titles and Labels
    title('SNR Distribution')
    xlabel('Signal to Noise Ratio')
    ylabel('Number of Observation')
    ax  = gca;
    ax.FontSize = 36;
    

% %     FileNString = strcat('Z:\SpecialNeeds\Sakib\Figures\', strcat(base(end-8:end-5),...
% %         base(end-3:end),'d', dates(date).name, 'band', num2str(band),'hist'));
% %     print(FileNString, '-djpeg','-r500')
%     hold off
%     else
%     end
%% Storing Data
     if Data_prc > 0.5
         SNR_95(date) = SNR_95prc;
         SNR_96(date) = SNR_96prc;
         SNR_97(date) = SNR_97prc;
         SNR_98(date) = SNR_98prc;
         SNR_99(date) = SNR_99prc;
         SNR_99_5(date) = SNR_99_5prc;
         SNR_99_7(date) = SNR_99_7prc;
         SNR_99_8(date) = SNR_99_8prc;
         SNR_99_9(date) = SNR_99_9prc;
     else
         SNR_95(date) = NaN;
         SNR_96(date) = NaN;
         SNR_97(date) = NaN;
         SNR_98(date) = NaN;
         SNR_99(date) = NaN;
         SNR_99_5(date) = NaN;
         SNR_99_7(date) = NaN;
         SNR_99_8(date) = NaN;
         SNR_99_9(date) = NaN;
     end

 end
     L8SNR.SNR_95{location, band} = SNR_95;
     L8SNR.SNR_96{location, band} = SNR_96;
     L8SNR.SNR_97{location, band} = SNR_97;
     L8SNR.SNR_98{location, band} = SNR_98;
     L8SNR.SNR_99{location, band} = SNR_99;
     L8SNR.SNR_99_5{location, band} = SNR_99_5;
     L8SNR.SNR_99_7{location,band} = SNR_99_7;
     L8SNR.SNR_99_8{location,band} = SNR_99_8;
     L8SNR.SNR_99_9{location,band} = SNR_99_9;
end
clearvars -except L8SNR
end
%% Storing all the data into matrix
% load('L8SNR_Blue_1prc.mat')

location = 17;
SNRall_96_mat = nan(location, 200);
SNRall_97_mat = nan(location, 200);
SNRall_98_mat = nan(location, 200);
SNRall_99_mat = nan(location, 200);
SNRall_99_5_mat = nan(location, 200);
SNRall_99_7_mat = nan(location, 200);
SNRall_99_8_mat = nan(location, 200);
SNRall_99_9_mat = nan(location, 200);

for location =  1:6 %17
    SNRall_96 = L8SNR.SNR_96{location, 1};
    SNRall_96_mat(location, 1:length(SNRall_96)) = SNRall_96;
    
    SNRall_97 = L8SNR.SNR_97{location, 1};
    SNRall_97_mat(location, 1:length(SNRall_97)) = SNRall_97;
    
    SNRall_98 = L8SNR.SNR_98{location, 1};
    SNRall_98_mat(location, 1:length(SNRall_98)) = SNRall_98;
    
    SNRall_99 = L8SNR.SNR_99{location, 1};
    SNRall_99_mat(location, 1:length(SNRall_99)) = SNRall_99;
    
    SNRall_99_5 = L8SNR.SNR_99_5{location, 1};
    SNRall_99_5_mat(location, 1:length(SNRall_99_5)) = SNRall_99_5;
    
    SNRall_99_7 = L8SNR.SNR_99_7{location, 1};
    SNRall_99_7_mat(location, 1:length(SNRall_99_7)) = SNRall_99_7;
    
    SNRall_99_8 = L8SNR.SNR_99_8{location, 1};
    SNRall_99_8_mat(location, 1:length(SNRall_99_8)) = SNRall_99_8;
    
    SNRall_99_9 = L8SNR.SNR_99_9{location, 1};
    SNRall_99_9_mat(location, 1:length(SNRall_99_9)) = SNRall_99_9;
end

%% Box Plot
band = 1;
L8_SNR_pre = {'370', '305', '230', '205'}; % At Typical Radiance Level
% L8_SNR_pre = {'1130', '1200', '950', '1000'}; % At High Radiance Level
band_name = {'BLUE', 'GREEN', 'RED', 'NIR'};
figure,
b = boxplot([SNRall_98_mat(:), SNRall_99_mat(:), SNRall_99_5_mat(:)],...
         'Labels', {'98%', '99%', '99.5%' },'Symbol','b.',  'OutlierSize', 30,'Whisker',1.5);
set(b,{'linew'},{2})
xLimits = xlim();
yLimits = ylim();
hold on
plot([xLimits(1) xLimits(2)], [str2double(L8_SNR_pre{band}) str2double(L8_SNR_pre{band})], 'g','LineWidth',3)
legend(strcat('L8 Pre-Launch SNR', {', '}, strcat(L8_SNR_pre{band})),'Location','northwest');
% ylim([0 1000])
xlabel('SNR at Percentage')
ylabel('Signal to Noise Ratio')
title('Signal to Noise Ratio of Landsat 8, Blue Band')
ax = gca;
ax.FontSize = 36;
grid on;
grid minor
ax.GridColor = 'k';
ax.MinorGridColor = 'k';

%% Individual value Plot 
    band = 4;
    figure(band)
    range = ones(1, length(SNRall_95_vec(:,3)));
    x = 0.2*range;
    plot(x, SNRall_95_vec(:,band), 'ro', 'MarkerSize', 10, 'LineWidth', 1);
    hold on;
    
    x = 0.4*range;
    plot(x, SNRall_96_vec(:,band), 'co', 'MarkerSize', 10, 'LineWidth', 1);
    hold on;
    
    x = 0.6*range;
    plot(x, SNRall_97_vec(:,band), 'bo', 'MarkerSize', 10, 'LineWidth', 1);
    hold on;
    
    x = 0.8*range;
    plot(x, SNRall_98_vec(:,band), 'go', 'MarkerSize', 10, 'LineWidth', 1);
    hold on; 
    
    x = 1.0*range;
    plot(x, SNRall_99_vec(:,band), 'mo', 'MarkerSize', 10, 'LineWidth', 1);
    
    hold on; 
    x = 1.2*range;
    plot(x, SNRall_99_5_vec(:,band), 'ko', 'MarkerSize', 10, 'LineWidth', 1);
    
    hold on; 
    x = 1.4*range;
    plot(x, SNRall_99_7_vec(:,band), 'mo', 'MarkerSize', 10, 'LineWidth', 1);
    
    % Set up axes.
    xlim([0, 1.4]);
    %ylim([0, 1000]);
    xlabel('SNR at Percentage')
    ylabel('Signal to Noise Ratio')
    title('Signal to Noise Ratio of L8')
    
    ax = gca;
    ax.XTick = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4];
    ax.XTickLabels = {'96%','97%','98%', '99%', '99.5%', '99.7%', '99.9%'};
    ax.FontSize = 36;
    grid on;
    grid minor
    ax.GridColor = 'k';
%   ax.GridAlpha = 0.8;
    ax.MinorGridColor = 'k';
    %ax.FontName = 'Times New Roman';
    
%% L8 SNR all band together
location_total = 17;
for band = 1%:4
       if band == 1
            load('L8SNR_Blue_2prc.mat')
            %load('L8SNR_125prc.mat')
            % Creating Nan matrix to store the values
            SNRall_96_mat = nan(location_total, 200); SNRall_97_mat = nan(location_total, 200);
            SNRall_98_mat = nan(location_total, 200); SNRall_99_mat = nan(location_total, 200);
            SNRall_99_5_mat = nan(location_total, 200); SNRall_99_7_mat = nan(location_total, 200);
            SNRall_99_8_mat = nan(location_total, 200); SNRall_99_9_mat = nan(location_total, 200);
            
        elseif band == 2
            load('L8SNR_Green_1prc.mat')
            % Creating Nan matrix to store the values
            SNRall_96_mat = nan(location_total, 200); SNRall_97_mat = nan(location_total, 200);
            SNRall_98_mat = nan(location_total, 200); SNRall_99_mat = nan(location_total, 200);
            SNRall_99_5_mat = nan(location_total, 200); SNRall_99_7_mat = nan(location_total, 200);
            SNRall_99_8_mat = nan(location_total, 200); SNRall_99_9_mat = nan(location_total, 200);   
            
        elseif band == 3
            load('L8SNR_Red_1prc.mat')
            % Creating Nan matrix to store the values
            SNRall_96_mat = nan(location_total, 200); SNRall_97_mat = nan(location_total, 200);
            SNRall_98_mat = nan(location_total, 200); SNRall_99_mat = nan(location_total, 200);
            SNRall_99_5_mat = nan(location_total, 200); SNRall_99_7_mat = nan(location_total, 200);
            SNRall_99_8_mat = nan(location_total, 200); SNRall_99_9_mat = nan(location_total, 200);
            
        elseif band == 4
            load('L8SNR_NIR_1prc.mat')
            % Creating Nan matrix to store the values
            SNRall_96_mat = nan(location_total, 200); SNRall_97_mat = nan(location_total, 200);
            SNRall_98_mat = nan(location_total, 200); SNRall_99_mat = nan(location_total, 200);
            SNRall_99_5_mat = nan(location_total, 200); SNRall_99_7_mat = nan(location_total, 200);
            SNRall_99_8_mat = nan(location_total, 200); SNRall_99_9_mat = nan(location_total, 200);
            
        end
    
    for location = 1:17
        SNRall_96 = L8SNR.SNR_96{location, 1};
        SNRall_96_mat(location, 1:length(SNRall_96)) = SNRall_96;
        
        SNRall_97 = L8SNR.SNR_97{location, 1};
        SNRall_97_mat(location, 1:length(SNRall_97)) = SNRall_97;
        
        SNRall_98 = L8SNR.SNR_98{location, 1};
        SNRall_98_mat(location, 1:length(SNRall_98)) = SNRall_98;

        SNRall_99 = L8SNR.SNR_99{location, 1};
        SNRall_99_mat(location, 1:length(SNRall_99)) = SNRall_99;

        SNRall_99_5 = L8SNR.SNR_99_5{location, 1};
        SNRall_99_5_mat(location, 1:length(SNRall_99_5)) = SNRall_99_5;
        
        SNRall_99_7 = L8SNR.SNR_99_7{location, 1};
        SNRall_99_7_mat(location, 1:length(SNRall_99_7)) = SNRall_99_7;
        
        SNRall_99_8 = L8SNR.SNR_99_8{location, 1};
        SNRall_99_8_mat(location, 1:length(SNRall_99_8)) = SNRall_99_8;
        
        SNRall_99_9 = L8SNR.SNR_99_9{location, 1};
        SNRall_99_9_mat(location, 1:length(SNRall_99_9)) = SNRall_99_9;
        
    end
        SNRall_96_vec(:, band) = SNRall_96_mat(:);
        SNRall_97_vec(:, band) = SNRall_97_mat(:);
        SNRall_98_vec(:, band) = SNRall_98_mat(:);
        SNRall_99_vec(:, band) = SNRall_99_mat(:);
        SNRall_99_5_vec(:, band) = SNRall_99_5_mat(:);
        SNRall_99_7_vec(:, band) = SNRall_99_7_mat(:);
        SNRall_99_8_vec(:, band) = SNRall_99_8_mat(:);
        SNRall_99_9_vec(:, band) = SNRall_99_9_mat(:);
 end

%% Ploting all band together
L8_SNR_pre = {'370', '305', '230', '205'}; % At Typical Radiance Level
%L8_SNR_pre = {'1130', '1200', '950', '1000'}; % At High Radiance Level
band_name = {'BLUE', 'GREEN', 'RED', 'NIR'};

for band = 1 %:3
    figure(band)
    if band == 1
    b = boxplot([SNRall_96_vec(:,band), SNRall_97_vec(:,band), SNRall_98_vec(:,band),...
        SNRall_99_vec(:,band), SNRall_99_5_vec(:,band)],...
         'Labels', {'96%','97%','98%','99%', '99.5%'},'Symbol','b.',  'OutlierSize',30,'Whisker',1.5);
    hold on
     
    elseif band == 2
        b = boxplot([SNRall_98_vec(:,band), SNRall_99_vec(:,band), SNRall_99_5_vec(:,band),...
            SNRall_99_7_vec(:,band), SNRall_99_8_vec(:,band)],...
             'Labels', {'98%','99%', '99.5%', '99.7%', '99.8%'},'Symbol','b.',  'OutlierSize',30,'Whisker',1.5);
         hold on
     
    elseif band == 3
        b = boxplot([SNRall_98_vec(:,band), SNRall_99_vec(:,band), SNRall_99_5_vec(:,band),...
            SNRall_99_7_vec(:,band), SNRall_99_8_vec(:,band)],...
             'Labels', {'98%','99%', '99.5%', '99.7%', '99.8%'},'Symbol','b.',  'OutlierSize',30,'Whisker',1.5);
         hold on
    end
    xLimits = xlim();
    yLimits = ylim();
    hold on
    plot([xLimits(1) xLimits(2)], [str2double(L8_SNR_pre{band}) str2double(L8_SNR_pre{band})], 'g','LineWidth',3)
    legend(strcat('L8 Pre-Launch SNR', {', '}, strcat(L8_SNR_pre{band})),'Location','northwest');
    hold on
    ylim([0 1000])
    %xlim([0 4.5])
    set(b,{'linew'},{2})
    xlabel('SNR at Percentage')
    ylabel('Signal to Noise Ratio')
    title(strcat('Signal to Noise Ratio of Landsat 8',  {', '}, strcat(band_name{band}, ' Band')))
    
    ax  = gca;
    grid on
    grid minor
    ax.FontSize = 36;
    hold off
end

%% Ploting Single SNR
    figure(1)
    b = boxplot([SNRall_98_vec(:,1), SNRall_99_5_vec(:,2), SNRall_99_5_vec(:,3)],...
         'Labels', {'Blue', 'Green', 'Red'},'Symbol','b.',  'OutlierSize',30,'Whisker',1.5);
    blue_q3 = quantile(SNRall_98_vec(:,1), 0.75);
    green_q3 = quantile(SNRall_99_vec(:,2), 0.75);
    hold on
    legend('L8 Pre-Launch SNR');
    ylim([0 1000])
    set(b,{'linew'},{2})
    %xlabel('BanSNR')
    ylabel('Signal to Noise Ratio')
    title('Signal to Noise Ratio of Landsat 8')
    ax  = gca;
    grid on
    grid minor
    ax.FontSize = 36;
    hold off
