function [Hyper_WLenNM, Hyper_Spectrum, Hy_WL_Sand,Hyper_Spectrum_all] = HyperData()
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
load('HyDessert_24Profile.mat');

% Storing both Vegitation and Sand Profile together 
DiffCol = size(Hyper_Spectrum,2) - size(HyDessert(2:end,:)',2);
Hy_WL_Sand = HyDessert(1,:)';

MeanRef2 = [HyDessert(2:end,:)' NaN(size(HyDessert(2:end,:)',1), DiffCol)];
Hyper_Spectrum_all = [Hyper_Spectrum; MeanRef2]; clear DiffCol MeanRef MeanRef2
end