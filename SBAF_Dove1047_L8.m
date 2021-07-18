clear all
%% SBAF- Dove 1047 to Landsat 8 with vegetation (SDSU Cal Site data) 
base_dir = 'Z:\SpecialNeeds\Sakib\SDSU_Cal_Site_Data\';
Directory = dir(base_dir); Directory([1 2]) = [];

set(0,'DefaultLegendAutoUpdate','off')
D1047_OLIRSR_plot(); hold on;

for file_no = 1: size(Directory)
    d = Directory(file_no).name;
    load([base_dir, d])
    
   if isfield(Out.Mod,'TOASun') == 0
       Hyper_Spectrum = NaN;
   else
       Hyper_Spectrum = Out.Mod.TotalRadS./Out.Mod.TOASun;
   end
    % Hyperspectral Wavelength
    Hyper_WLenNM = Out.Mod.WaveS;
    Hyper_WLenNM_range = Hyper_WLenNM;
    Hyper_WLenNM_range(Hyper_WLenNM_range >1000) = nan;
   % Hyper_Spectrum(file_no,:) = Hyper_Spectrum;
    plot(Hyper_WLenNM_range, Hyper_Spectrum)
    hold on
    %%% Banding Reference Sensor OLI- Satnum = 24
   % [labels_L8, HyBanded_L8] = bander(Hyper_WLenNM, Hyper_Spectrum, 24);

    %%% Banding Target Sensor Dove R(1047)- Satnum = 84
   % [labels_1047, HyBanded_1047] = bander(Hyper_WLenNM, Hyper_Spectrum, 84);
   % SBAF_L8_D1047_SCalSite(file_no, :) = HyBanded_L8(:, 2:5)./HyBanded_1047;
end
%xlim([400 1000])
%ylim([0 0.002])

%% Different Types of Sand Profile
load('HyDessert_24Profile.mat');
set(0,'DefaultLegendAutoUpdate','off')
D1047_OLIRSR_plot(); hold on;
plot(HyDessert(1,:), HyDessert(2:end,:))

for i = 1:size(HyDessert,1)-1
    Hyper_Spectrum = HyDessert(i+1,:);
   %%% Banding Reference Sensor OLI- Satnum = 24
    [labels_L8, HyBanded_L8] = bander(HyDessert(1,:), Hyper_Spectrum, 24);

   %%% Banding Target Sensor Dove R(1047)- Satnum = 84
    [labels_1047, HyBanded_1047] = bander(HyDessert(1,:), Hyper_Spectrum, 84);
    SBAF_L8_D1047_SCalSite(i, :) = HyBanded_L8(:, 2:5)./HyBanded_1047;
end 

%% Sand PICS profile-- 
%%% Loading PICS Hyperspectral Data
load('hyperion_over_Libya_1_view_solar.mat'); meanRef_L1 = meanReflectance;
load('hyperion_over_Niger_1_view_solar.mat'); meanRef_N1 = meanReflectance;
%load('hyperion_over_Niger_1_view_solar.mat'); meanRef_N2 = meanReflectance;
load('hyperion_over_Sudan_1_view_solar.mat'); meanRef_S1 = meanReflectance; 
MeanRef = [meanRef_L1; meanRef_N1(2:end,:); meanRef_S1(2:end,:)];
clear meanRef_L1 meanRef_N1 meanRef_S1

%%% Plotting Dove 1047 and OLI RSR with Sand Hyper Profile
set(0,'DefaultLegendAutoUpdate','off')
D1047_OLIRSR_plot(); hold on;
plot(MeanRef(1,:), MeanRef(2:end,:))

%% Calculating SBAF for Sand: OLI and Dove 1047
clear all
%% Loading the Hyperspectral profiles and OLI Reflectance
[Hyper_WLenNM, Hyper_Spectrum, Hy_WL_Sand, Hyper_Spectrum_all] = HyperData();
L8_allMean = OLI_DoveTOACCdata(); 
ref = 1;
OLI_ref = L8_allMean(ref,:); 
[labels_L8, HyBanded_L8] = bander(Hyper_WLenNM, Hyper_Spectrum(:,1), 24);

%% Plotting the Original
figure, plot(OLI_ref,'r.', 'MarkerSize', 50); hold on
plot(HyBanded_L8(2:5), 'bp', 'LineWidth',2,'MarkerSize', 30)
ylim([0 0.7])
x = {'Blue', '', 'Green', '', 'Red', '', 'NIR'};
set(gca,'xticklabel', x)
title('Landsat 8 and Hyperspectral Banded Reflectance');
ylabel('TOA Reflectance')
xlabel('Bands')
hold on; grid on; grid minor; ax = gca; ax.FontSize = 30; ax.GridColor = 'k';
ax.MinorGridColor = 'k'; % ax.GridAlpha = 0.8;
legend({'Landsat 8', 'Hyperspectral Banded'}, 'Location','northwest')

%% Scaling and ploting TOA Reflectance
M_vegi = mean(OLI_ref./HyBanded_L8(2:5));
figure, plot(OLI_ref,'r.', 'MarkerSize', 50); hold on
plot(M_vegi.*HyBanded_L8(2:5), 'bp', 'LineWidth',2,'MarkerSize', 30)
ylim([0 0.7])
x = {'Blue', '', 'Green', '', 'Red', '', 'NIR'};
set(gca,'xticklabel', x)
title('Landsat 8 and Hyperspectral Banded and Scaled Reflectance');
ylabel('TOA Reflectance')
xlabel('Bands')
hold on; grid on; grid minor; ax = gca; ax.FontSize = 30; ax.GridColor = 'k';
ax.MinorGridColor = 'k'; % ax.GridAlpha = 0.8;
legend({'Landsat 8', 'Hyperspectral Banded'}, 'Location','northwest')
%%
% figure, plot(OLI_ref, 'r.', 'MarkerSize', 40); hold on
% plot(HyB_L8_m, 'bp', 'LineWidth',2,'MarkerSize', 20)
% figure, plot(OLI_ref,'.'); hold on
% plot(M_vegi.*HyBanded_L8(2:5), '*')
x = {'Blue', '', 'Green', '', 'Red', '', 'NIR'};
set(gca,'xticklabel', x)
title('Landsat 8 and Hyperspectral Banded Reflectance');
ylabel('TOA Reflectance')
xlabel('Bands')
hold on; grid on; grid minor; ax = gca; ax.FontSize = 30; ax.GridColor = 'k';
ax.MinorGridColor = 'k'; % ax.GridAlpha = 0.8;
legend({'Landsat 8', 'Hyperspectral Banded'}, 'Location','northwest')

%% Plotting SBAF Dove 1047 and Landsat 8
figure,
for p = 1: size(SBAF_L8_D1047_SCalSite, 1)
     plot(SBAF_L8_D1047_SCalSite(p,:), 'g.', 'MarkerSize', 40); %ylim([0.97 1.03]);
     hold on
end

%ylim([0.94 1.04])
x = {'Blue', '', 'Green', '', 'Red', '', 'NIR'};
set(gca,'xticklabel', x)

title('SBAF for Dove(1047) to Landsat 8 over Desert');
ylabel('SBAF_{D1047-to-Landsat8}')
xlabel('Bands')
hold on; grid on; grid minor; ax = gca; ax.FontSize = 30;ax.GridColor = 'k';
ax.MinorGridColor = 'k'; % ax.GridAlpha = 0.8;



