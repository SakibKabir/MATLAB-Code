clear all
%% Dove Data loading and organizing for plot
%clear all
%load('DoveSNR_125prc.mat')
location_total = 17;
% Creating Nan matrix to store the values
SNRall_95_mat = nan(location_total, 200);
SNRall_96_mat = nan(location_total, 200); SNRall_97_mat = nan(location_total, 200);
SNRall_98_mat = nan(location_total, 200); SNRall_99_mat = nan(location_total, 200);
SNRall_99_5_mat = nan(location_total, 200); SNRall_99_7_mat = nan(location_total, 200);
SNRall_99_8_mat = nan(location_total, 200); SNRall_99_9_mat = nan(location_total, 200);
   
for band = 1:4
   
for location = 1:17
        SNRall_95 = DoveSNR.SNR_95_D1047{location, band};
        SNRall_95_mat(location, 1:length(SNRall_95)) = SNRall_95;
        
        SNRall_96 = DoveSNR.SNR_96_D1047{location, band};
        SNRall_96_mat(location, 1:length(SNRall_96)) = SNRall_96;
    
        SNRall_97 = DoveSNR.SNR_97_D1047{location, band};
        SNRall_97_mat(location, 1:length(SNRall_97)) = SNRall_97;
        
        SNRall_98 = DoveSNR.SNR_98_D1047{location, band};
        SNRall_98_mat(location, 1:length(SNRall_98)) = SNRall_98;

        SNRall_99 = DoveSNR.SNR_99_D1047{location, band};
        SNRall_99_mat(location, 1:length(SNRall_99)) = SNRall_99;

        SNRall_99_5 = DoveSNR.SNR_99_5_D1047{location, band};
        SNRall_99_5_mat(location, 1:length(SNRall_99_5)) = SNRall_99_5;
        
        SNRall_99_7 = DoveSNR.SNR_99_7_D1047{location, band};
        SNRall_99_7_mat(location, 1:length(SNRall_99_7)) = SNRall_99_7;
        
        SNRall_99_8 = DoveSNR.SNR_99_8_D1047{location, band};
        SNRall_99_8_mat(location, 1:length(SNRall_99_8)) = SNRall_99_8;
        
        SNRall_99_9 = DoveSNR.SNR_99_9_D1047{location, band};
        SNRall_99_9_mat(location, 1:length(SNRall_99_9)) = SNRall_99_9;
end

    SNRall_95_vec(:, band) = SNRall_95_mat(:);
    SNRall_96_vec(:, band) = SNRall_96_mat(:);
    SNRall_97_vec(:, band) = SNRall_97_mat(:);
    SNRall_98_vec(:, band) = SNRall_98_mat(:);
    SNRall_99_vec(:, band) = SNRall_99_mat(:);
    SNRall_99_5_vec(:, band) = SNRall_99_5_mat(:);
    SNRall_99_7_vec(:, band) = SNRall_99_7_mat(:);
    SNRall_99_8_vec(:, band) = SNRall_99_8_mat(:);
    SNRall_99_9_vec(:, band) = SNRall_99_9_mat(:);
end 
%% Ploting all Dove band together
band_name = {'BLUE', 'GREEN', 'RED', 'NIR'};

for band = 3 %:3
    figure(band)
    if band == 1
    b = boxplot([SNRall_97_vec(:,band), SNRall_98_vec(:,band),...
        SNRall_99_vec(:,band), SNRall_99_5_vec(:,band)],...
         'Labels', {'97%','98%','99%', '99.5%'},'Symbol','b.',  'OutlierSize',30,'Whisker',1.5);
    hold on
     
    elseif band == 2
        b = boxplot([SNRall_97_vec(:,band),SNRall_98_vec(:,band),...
            SNRall_99_vec(:,band), SNRall_99_5_vec(:,band)],...
             'Labels', {'97%','98%','99%', '99.5%'},'Symbol','b.',  'OutlierSize',30,'Whisker',1.5);
         hold on
     
    elseif band == 3
        b = boxplot( [SNRall_97_vec(:,band),SNRall_98_vec(:,band),...
            SNRall_99_vec(:,band), SNRall_99_5_vec(:,band)],...
             'Labels', { '97%','98%', '99%', '99.5%'},'Symbol','b.',  'OutlierSize',30,'Whisker',1.5);
         hold on
         
    elseif band == 4
        b = boxplot([SNRall_97_vec(:,band),SNRall_98_vec(:,band),...
            SNRall_99_vec(:,band), SNRall_99_5_vec(:,band)],...
             'Labels', { '98%', '99%', '99.5%', '99.7%', '99.8%', '99.9%'},'Symbol','b.',  'OutlierSize',30,'Whisker',1.5);
         hold on
    end
    xLimits = xlim();
    yLimits = ylim();
    hold on
    %plot([xLimits(1) xLimits(2)], [str2double(L8_SNR_pre{band}) str2double(L8_SNR_pre{band})], 'g','LineWidth',3)
    %legend(strcat('L8 Pre-Launch SNR', {', '}, strcat(L8_SNR_pre{band})),'Location','northwest');
   % hold on
    ylim([0 100])
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

%% Landsat 8 Data loading and organizing for plot
%clear all
%load('L8SNR_125prc.mat')
load('L8SNR_250prc.mat')
location_total = 17;
% Creating Nan matrix to store the values
SNRall_95_mat = nan(location_total, 200);
SNRall_96_mat = nan(location_total, 200); SNRall_97_mat = nan(location_total, 200);
SNRall_98_mat = nan(location_total, 200); SNRall_99_mat = nan(location_total, 200);
SNRall_99_5_mat = nan(location_total, 200); SNRall_99_7_mat = nan(location_total, 200);
SNRall_99_8_mat = nan(location_total, 200); SNRall_99_9_mat = nan(location_total, 200);
   
for band = 1:4
    L8_band = band + 1;
for location = 1:17
        SNRall_95 = L8SNR.SNR_95{location, L8_band};
        SNRall_95_mat(location, 1:length(SNRall_95)) = SNRall_95;
        
        SNRall_96 = L8SNR.SNR_96{location, L8_band};
        SNRall_96_mat(location, 1:length(SNRall_96)) = SNRall_96;
    
        SNRall_97 = L8SNR.SNR_97{location, L8_band};
        SNRall_97_mat(location, 1:length(SNRall_97)) = SNRall_97;
        
        SNRall_98 = L8SNR.SNR_98{location, L8_band};
        SNRall_98_mat(location, 1:length(SNRall_98)) = SNRall_98;

        SNRall_99 = L8SNR.SNR_99{location, L8_band};
        SNRall_99_mat(location, 1:length(SNRall_99)) = SNRall_99;

        SNRall_99_5 = L8SNR.SNR_99_5{location, L8_band};
        SNRall_99_5_mat(location, 1:length(SNRall_99_5)) = SNRall_99_5;
        
        SNRall_99_7 = L8SNR.SNR_99_7{location, L8_band};
        SNRall_99_7_mat(location, 1:length(SNRall_99_7)) = SNRall_99_7;
        
        SNRall_99_8 = L8SNR.SNR_99_8{location, L8_band};
        SNRall_99_8_mat(location, 1:length(SNRall_99_8)) = SNRall_99_8;
        
        SNRall_99_9 = L8SNR.SNR_99_9{location, L8_band};
        SNRall_99_9_mat(location, 1:length(SNRall_99_9)) = SNRall_99_9;
end
    SNRall_95_vec(:, band) = SNRall_95_mat(:);
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

% Landsat bands starts from Band 1 = blue
for band = 1 :4
    figure(band)
    if band == 1
    b = boxplot([SNRall_96_vec(:,band), SNRall_97_vec(:,band), SNRall_98_vec(:,band),...
        SNRall_99_vec(:,band), SNRall_99_5_vec(:,band)],...
         'Labels', {'96%','97%','98%','99%', '99.5%'},'Symbol','b.',  'OutlierSize',30,'Whisker',1.5);
    hold on
     
    elseif band == 2
        b = boxplot([SNRall_96_vec(:,band), SNRall_97_vec(:,band),SNRall_98_vec(:,band),...
            SNRall_99_vec(:,band), SNRall_99_5_vec(:,band)],...
             'Labels', {'96%','97%','98%','99%', '99.5%'},'Symbol','b.',  'OutlierSize',30,'Whisker',1.5);
         hold on
     
    elseif band == 3
        b = boxplot([ SNRall_98_vec(:,band),...
            SNRall_99_vec(:,band), SNRall_99_5_vec(:,band), SNRall_99_7_vec(:,band), SNRall_99_8_vec(:,band)],...
             'Labels', { '98%', '99%', '99.5%', '99.7%', '99.8%'},'Symbol','b.',  'OutlierSize',30,'Whisker',1.5);
         hold on
         
    elseif band == 4
        b = boxplot([ SNRall_98_vec(:,band),...
            SNRall_99_vec(:,band), SNRall_99_5_vec(:,band), SNRall_99_7_vec(:,band), SNRall_99_8_vec(:,band), SNRall_99_9_vec(:,band)],...
             'Labels', { '98%', '99%', '99.5%', '99.7%', '99.8%', '99.9%'},'Symbol','b.',  'OutlierSize',30,'Whisker',1.5);
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

%% Ploting Single SNR all bands
L8_SNR_pre = {'370', '305', '230', '205'}; % At Typical Radiance Level
band_colors={'c','g','r','m'};
band_name = {'BLUE', 'GREEN', 'RED', 'NIR'};

figure(1)
b = boxplot([SNRall_97_vec(:,1), SNRall_99_5_vec(:,2), SNRall_99_7_vec(:,3)],...
     'Labels', {'Blue', 'Green', 'Red'},'Symbol','b.',  'OutlierSize',30,'Whisker',1.5);
hold on
for band = 1:3
  plot(band, str2double(L8_SNR_pre{band}),'.' ,'color', band_colors{band},'markers', 36 )
  hold on
end

legend(strcat('L8 Pre-Launch SNR', {', '}, strcat(L8_SNR_pre{band})),'Location','northeast');
hold on

set(gca,'box','off')
h = copyobj(gca,gcf);
plot(h, 1:3,[370, 305, 230])

xLimits = xlim();
yLimits = ylim();
    
for band = 1:4
    plot([xLimits(1) xLimits(2)], [str2double(L8_SNR_pre{band}) str2double(L8_SNR_pre{band})], 'color', band_colors{band},'LineWidth',3)
    legend(strcat('L8 Pre-Launch SNR', {', '}, strcat(L8_SNR_pre{band})),'Location','northwest');
    hold on
end

% legend('L8 Pre-Launch SNR');
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
