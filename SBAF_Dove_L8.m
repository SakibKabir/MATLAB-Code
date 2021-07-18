clear all
%%% Learning How to perform SBAF- Ref OLI and Target MSI S2B
load('Z:\SpecialNeeds\Sakib\From Bipin\hyperion_over_sudan_1.mat')

Hy_WLenNM = meanReflectance(1,:);
Hy_Spectrum = meanReflectance(2:end,:);

%%% Reference Sensor- OLI (Satnum 24)
for sel = 1:168
    [labels, Banded] = bander(Hy_WLenNM, Hy_Spectrum(sel,:), 24);
    bandedHyRef_L8(sel,:) = Banded;
end
%figure, plot(bandedHyRef_L8(:,1), '*')

%% Target Sensor- S2B MSI (Satnum 36)
for sel = 1:168
    [labels, Banded] = bander(Hy_WLenNM, Hy_Spectrum(sel,:), 36);
    bandedHyRef_S2B(sel,:) = Banded;
end

plot(bandedHyRef_S2B(:,1), '.')

for band = 1:4
   SBAF_bands(band) = bandedHyRef_L8(1, band)/bandedHyRef_S2B(1, band);
end

band = 5;
S_band = 13;
SBAF_bands(band) = bandedHyRef_L8(1, band)/bandedHyRef_S2B(1, S_band);

figure, plot(SBAF_bands, 'h', 'MarkerSize', 10)

%% Checking the SBAF of Sonoron Desert from Morakots file
% Sonoron Desert- path 38 and row 38.
%%% Sonoran Hyperion profile from Morakot
load('hyperion profile for sonoran.mat')
Hy_Spectrum_SoD = Sonoranhyperionprofile; 

%%% Sudan 1 Hyperion for Wavelengths
load('Z:\SpecialNeeds\Sakib\From Bipin\hyperion_over_sudan_1.mat')
Hy_WLenNM = meanReflectance(1,:);

%%% Banding Reference Sensor OLI- Satnum = 24
[labels_L8, HyBanded_L8] = bander(Hy_WLenNM, Hy_Spectrum_SoD, 24);

%%% Checking the calculation with S2B
[labels, Banded_S2B] = bander(Hy_WLenNM, Hy_Spectrum_SoD, 36);
SBAF_L8_S2B = HyBanded_L8(:,1:4)./Banded_S2B(:,1:4);
%figure, plot(SBAF_L8_S2B, '.', 'MarkerSize', 40); ylim([0.96 1.03]);

%%% Banding Target Sensor Dove R(1058)- Satnum = 83
[labels_1058, HyBanded_1058] = bander(Hy_WLenNM, Hy_Spectrum_SoD, 83);
SBAF_L8_D1058 = HyBanded_L8(:, 2:5)./HyBanded_1058;
figure, plot(SBAF_L8_D1058, 'h','MarkerFaceColor','b', 'MarkerSize', 20); ylim([0.96 1.03]);

%% Banding Target Sensor Dove R(1047)- Satnum = 84
[labels_1047, HyBanded_1047] = bander(Hy_WLenNM, Hy_Spectrum_SoD, 84);
SBAF_L8_D1047 = HyBanded_L8(:, 2:5)./HyBanded_1047;
figure, plot(SBAF_L8_D1047, 'h','MarkerFaceColor','b', 'MarkerSize', 20); ylim([0.96 1.03]);

%% Printing L8 and Dove-R 1058 RSR with Hyperion
File_L8 = 'Z:\SpecialNeeds\Sakib\Miscellaneous\Ball_BA_RSR.xlsx';
Dfile_1058 = 'Z:\ImageDrive\PlanetLabs\Processed\1047\Preview\RSRs\1047.csv';
RSR_1058 = csvread(Dfile_1058, 1, 0);
% legend 
lege_1058_L8 = {'1058 Blue','L8 Blue','1058 Green', 'L8 Green','1058 Red','L8 Red', '1058 NIR', 'L8 NIR', 'S.Desert Hyperion'};

%%% Bands and Colors
bands ={'Blue', 'Green', 'Red', 'NIR'};
band_colors={'b','g','r','m'};

for band =1: 4
%%% Dove 1058 RSR
    plot(RSR_1058(:,1), RSR_1058(:,band+1), 'color', band_colors{band}, 'LineStyle', '-', 'LineWidth', 2.5)
    hold on
    
%%% Landsat 8 RSR 
    L8_RSR = xlsread(File_L8, band+1);
    plot(L8_RSR(:,1), L8_RSR(:,2), 'color', band_colors{band}, 'LineStyle', '--','LineWidth', 1)
    hold on
end

hold on
plot(Hy_WLenNM, Hy_Spectrum_SoD, 'LineWidth', 1.5 ); 
xlim([400 1000]);
ylim([0 1.1]);

legend(lege_1058_L8,'FontSize',18, 'Location','northeast');

title('Dove-R(1058) and L8 Relative Spectral Response with S.Desert Hyperion');
ylabel('Relative Response')
xlabel('Wavelength (nm)')

hold on
grid on
grid minor
ax  = gca;
ax.FontSize = 30;
ax.GridColor = 'k';
%   ax.GridAlpha = 0.8;
ax.MinorGridColor = 'k';

x={'Blue', '', 'Green', '', 'Red', '', 'NIR'};
set(gca,'xticklabel', x)

title('SBAF for Dove-R(1058) to Landsat 8');
ylabel('SBAF_{DR1058toLandsat8}')
xlabel('Bands')

%%
base='Z:\ImageDrive\Hyperion\EO1\P038\R038';
dates = dir(base); dates([1 2])=[];
date = '20170124';
Directory = dir(fullfile (base, date));
Directory([1 2]) = [];

%MTL File Parsing for Landsat 7
MTL=dir(fullfile(base, date, Directory.name,'*MTL.txt'));
[MTL_List_L8, vaule]= MTL_parser_L8(fullfile(MTL.folder, MTL.name));

% Dove File 
Dove_tif_name = '20190218_180340_16_1058_3B_AnalyticMS.tif';

%Extracting all the values from Metadata file
MData_file=fullfile(base, date, Directory.name, '20190218_180340_16_1058_3B_AnalyticMS_metadata.xml');
[MData_values]= xml2struct_new_v(MData_file);
[MData_values_old]= xml2struct(MData_file);

