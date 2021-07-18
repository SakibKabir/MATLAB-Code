clear all

if ispc 
    base_dir = 'Z:\SpecialNeeds\Sakib\SDSU_Cal_Site_Data\';
elseif isunix
    base_dir = '/home/sakib.kabir/zdrive/SpecialNeeds/Sakib/SDSU_Cal_Site_Data/';
end

Directory = dir(base_dir); Directory([1 2]) = [];

%%% Vegitation 
for file_no = 1: size(Directory)
    d = Directory(file_no).name;
    load([base_dir, d]); clear d
    
   if isfield(Out.Mod,'TOASun') == 0
       Hyper_Spectrum(:,file_no) = NaN; 
   else
       Hyper_Spectrum(:,file_no) = Out.Mod.TotalRadS./Out.Mod.TOASun;
   end
       Hyper_Spectrum(Hyper_Spectrum <=0) =nan;
       Hyper_Spectrum(Hyper_Spectrum == Inf) =nan;
      % Hyperspectral Wavelength
      Hyper_WLenNM(:,1) = Out.Mod.WaveS; clear Out
end

%% Sand
load('hyperion_over_Libya_1_view_solar.mat'); meanRef_L1 = meanReflectance;
load('hyperion_over_Niger_1_view_solar.mat'); meanRef_N1 = meanReflectance;
load('hyperion_over_Niger_1_view_solar.mat'); meanRef_N2 = meanReflectance;
load('hyperion_over_Sudan_1_view_solar.mat'); meanRef_S1 = meanReflectance; 
MeanRef = [meanRef_L1; meanRef_N1(2:end,:); meanRef_S1(2:end,:)];
clearvars -except Hyper_Spectrum Hyper_WLenNM MeanRef

% Storing both Vegitation and Sand Profile together 
DiffCol = size(Hyper_Spectrum,2) - size(MeanRef(2:end,:)',2);
Hy_WL_Sand = MeanRef(1,:)';

MeanRef2 = [MeanRef(2:end,:)' NaN(size(MeanRef(2:end,:)',1), DiffCol)];
Hyper_Spectrum_all = [Hyper_Spectrum; MeanRef2]; clear DiffCol MeanRef MeanRef2

%% Calculating SBAF with Appropiate Profile
%%% Only need the OLI Mean Reflectance
L8_allMean = OLI_DoveTOACCdata();

for ref = 1: size(L8_allMean) % all the reflectances inside the ROI
    for NoPro = 1:size(Hyper_Spectrum_all, 2)
        
        %%% For Vegitation
        Hy_Spectrum_Vegi = Hyper_Spectrum_all(1:size(Hyper_Spectrum,1), NoPro);
        %Hy_WL_Vegi = Hyper_WLenNM;
        %%% Banding Reference Sensor: OLI- Satnum = 24
        [Labels_L8_Vegi, HyBanded_L8_Vegi] = bander(Hyper_WLenNM, Hy_Spectrum_Vegi, 24);
        SSE(1, NoPro) = sse(L8_allMean(ref,:), HyBanded_L8_Vegi(2:5)); 
        clear Hy_Spectrum_Vegi Labels_L8_Vegi HyBanded_L8_Vegi 
        
        %%% For Sand
        Hy_Spectrum_Sand = Hyper_Spectrum_all(size(Hyper_Spectrum,1)+1:end, NoPro);
        %%% Banding Reference Sensor: OLI- Satnum = 24
        [Labels_L8_Sand, HyBanded_L8_sand] = bander(Hy_WL_Sand, Hy_Spectrum_Sand, 24);
        SSE(2, NoPro) = sse(L8_allMean(ref,:), HyBanded_L8_sand(2:5));  
        clear Hy_Spectrum_Sand Labels_L8_Sand HyBanded_L8_sand 
    end; clear NoPro       
    
 SSE(SSE==0) = nan; [r, c]= find(SSE == min(min(SSE))); % clear SSE
%%
if r(1) == 1 % Vegi
    Hyper = Hyper_Spectrum_all(1:size(Hyper_Spectrum, 1), c(1));
    Hyper_WL = Hyper_WLenNM;

elseif r(1) == 2 % Sand
    Hyper = Hyper_Spectrum_all(size(Hyper_Spectrum, 1)+1:end, c(1));
    Hyper_WL = Hy_WL_Sand;
end 

    %%% Banding Reference Sensor OLI- Satnum = 24
    [labels_L8, HyBanded_L8] = bander(Hyper_WL, Hyper, 24);
    %%% Banding Target Sensor Dove R(1047)- Satnum = 84
    [labels_1047, HyBanded_1047] = bander(Hyper_WL, Hyper, 84);
    SBAF_L8_D1047_SCalSite(ref, :) = HyBanded_L8(:, 2:5)./HyBanded_1047;
end 

clearvars -except SBAF_L8_D1047_SCalSite
save('SBAF_L8_D1047_SCalSite.mat')

%% Plotting the difference in OLI and Hyper banded Reflectance
% figure, plot(L8_allMean(ref,:), '.', 'MarkerSize', 50)
% hold on; plot(HyBanded_L8_Vegi(2:5), 'h','MarkerFaceColor', 'r', 'MarkerSize', 25); ylim([0 0.5])
% %hold on; plot(HyBanded_L8(2:5), 'h','MarkerFaceColor', 'r', 'MarkerSize', 25); ylim([0 0.5])
% x = {'Blue', '', 'Green', '', 'Red', '', 'NIR'}; set(gca,'xticklabel', x)
% title('OLI and Hyperspectral Banded Refelctance');
% xlabel('Bands'); ylabel('TOA Reflectance')
% legend('Landsat 8', 'Hyperspectral Banded')
% hold on; grid on; grid minor; ax  = gca; ax.FontSize = 30; ax.GridColor = 'k';