
bands=  {'1' '2' '3' '4' '5' '6' '7' '8'}; 
% %path/row- 
base='Z:\SpecialNeeds\Sakib\DoveData\April Data\Dove-R Data\1058\20190218\P038\R038';
dates = dir(base); dates([1 2])=[];
date = '20190218';
Directory = dir(fullfile (base, date));
Directory([1 2]) = [];

%MTL File Parsing for Landsat 7
MTL=dir(fullfile(base, date, Directory.name,'*MTL.txt'));
[MTL_List_L8, vaule]= MTL_parser_L8(fullfile(MTL.folder, MTL.name));

% Dove File Path
DoveDataPath = fullfile(base, date, Directory.name, '*_AnalyticMS.tif');
DoveDirectory = dir(DoveDataPath);
DoveDataPathXML = fullfile(base, date, Directory.name, '*_metadata.xml');
DoveDirectoryXML = dir(DoveDataPathXML);

%%% Path/Row - 038/037- 3 Dove Scene
RefMat_dove_16 = [0 -3; 3 0; 701148 3642078]; % for scene 16
RefMat_dove_22 = [0 -3; 3 0; 698382 3627618]; % for scene 16
RefMat_dove_29 = [0 -3; 3 0; 695658 3613158]; % for scene 29
% RefMat_dove = {RefMat_dove_16, RefMat_dove_22, RefMat_dove_29};

%%% Path/Row - 038/038 - 4 Dove Scene
RefMat_dove_36 = [0 -3; 3 0; 692901 3598653]; % for scene 36
RefMat_dove_42 = [0 -3; 3 0; 690126 3584151]; % for scene 42
RefMat_dove_49 = [0 -3; 3 0; 687342 3569691]; % for scene 49
RefMat_dove_55 = [0 -3; 3 0; 684600 3555240]; % for scene 55
RefMat_dove = {RefMat_dove_36, RefMat_dove_42, RefMat_dove_49, RefMat_dove_55};

% Dove Scenes
 for scene = 1:4

    %Extracting all the values from Metadata file
    MData_file = fullfile(base, date, Directory.name, DoveDirectoryXML(scene).name);
    [MData_values]= xml2struct_new_v(MData_file);
    [MData_values_old]= xml2struct(MData_file);
    
for D_band =1%:4
    L8_band = D_band+1;
    L8_Folder_info=dir(fullfile(base, date, Directory.name, strcat(Directory.name,'_B', bands{L8_band},'.tif')));  
    L8_image_file= fullfile(L8_Folder_info.folder, L8_Folder_info.name);
    L8_image_info = geotiffinfo(L8_image_file);

    [DN_L8, R_L8] = geotiffread(L8_image_file);
    DN_L8 = double(DN_L8);
    DN_L8(DN_L8 == 0)= nan;

    % Mulitplicative and additive factors
    rmb= MTL_List_L8.RADIOMETRIC_RESCALING.(strcat('RADIANCE_MULT_BAND_', bands{L8_band}));
    rab= MTL_List_L8.RADIOMETRIC_RESCALING.(strcat('RADIANCE_ADD_BAND_', bands{L8_band}));
    % L8 TOA Radiance
    L8_TOArad = DN_L8*rmb+rab;

    %%% Dove
    % Dove Image file- scene by scene
    Dove_image_file = fullfile(base, date, Directory.name, DoveDirectory(scene).name);

    % Reading Dove Image
    Dove_image_all_band =imread(Dove_image_file);
    DN_Dove = Dove_image_all_band(:,:, D_band);
    DN_Dove = double(DN_Dove);
    DN_Dove(DN_Dove == 0) = NaN;
    band_sf=0.01;
    Dove_TOArad = DN_Dove*band_sf;
end

%% Dove
% Checking Dove Radiance Image
% figure, imagesc(Dove_TOArad); colorbar
Dove_TOArad_resampled = imresize(Dove_TOArad, 3/30); % downsampling
% figure, imagesc(Dove_TOArad_resampled); colorbar
% New corner Coordinates

[nR,nC] = map2pix(R_L8, RefMat_dove{scene}(3,1), RefMat_dove{scene}(3,2));
nR=round(nR);
nC=round(nC);
[sR, sC]= size(Dove_TOArad_resampled); 
% Creating NaN image same size as Landsat 8 Image
D_TOArad_final = NaN(size(L8_TOArad,1), size(L8_TOArad, 2)); %size is same for all L8 band
% Inserting the Dove resampled dove image to NaN image
D_TOArad_final(nR:nR+sR-1,nC:nC+sC-1)= Dove_TOArad_resampled;
% figure, imagesc(D_TOArad_final); colorbar

D_TOArad_final_filter = D_TOArad_final;
D_TOArad_final_filter(~isnan(D_TOArad_final_filter)) = 1;
% figure, imagesc(D_TOArad_final_filter); colorbar

%% Landsat 8
%%% Checking L8 Radiance Image
%figure, imagesc(L8_TOArad); colorbar
L8_TOArad_DROI = L8_TOArad./ D_TOArad_final_filter;
% figure, imagesc(L8_TOArad_DROI); colorbar
L8_TOArad_filtered = L8_TOArad_DROI; % for filtering TOA Rad with typical
L8_TOArad_filtered(L8_TOArad_filtered < 39.5 | L8_TOArad_filtered > 40) = NaN;
% figure, imagesc(L8_TOArad_filtered); colorbar
% figure, histogram(L8_TOArad_filtered) 

%%% Checking L8 Standard Deviation image
window = 3;
%Dove_mean_image = conv2(Dove_TOArad, ones(window)./(window*window),'same');
k= (window-1)/2; 
%k_L8= (window_L8-1)/2; 
L8_SD_image = movingstd2(L8_TOArad_DROI, k);
% figure, imagesc(L8_SD_image); colorbar
L8_SD_image_filtered = L8_SD_image;
L8_SD_image_filtered(L8_SD_image_filtered < 0.10 | L8_SD_image_filtered> 0.12) = NaN; % 0.11+/-0.01
% figure, imagesc(L8_SD_image_filtered); colorbar
% figure, histogram(L8_SD_image_filtered) 

L8_SNR_image_filtered = L8_TOArad_filtered./L8_SD_image_filtered;
% figure, imagesc(L8_SNR_image_filtered); colorbar
% figure, histogram(L8_SNR_image_filtered) 
% figure, histogram(L8_SNR_image_filtered, 100) 
L8_max_SNR = max(max(L8_SNR_image_filtered));
L8_min_SNR = min(min(L8_SNR_image_filtered));
% sum(sum(~isnan(L8_SNR_image_filtered)))

%%% Create a Mask for dove data filtering with 1 and NaN
% With Rad value
L8_Rad_image_mask = L8_TOArad_filtered;
L8_Rad_image_mask(~(isnan(L8_Rad_image_mask))) = 1;
%figure, imagesc(L8_Rad_image_mask); colorbar

% With SNR value
L8_SNR_image_mask = L8_SNR_image_filtered;
L8_SNR_image_mask(~(isnan(L8_SNR_image_mask))) = 1;
%figure, imagesc(L8_SNR_image_mask); colorbar

%% ===================================================

%%%% This Section is for both Dove and L8
k= (window-1)/2; 
%k_L8= (window_L8-1)/2; 
D_SD_image = movingstd2(D_TOArad_final, k);
% figure, imagesc(D_SD_image); colorbar

% Dove Rad pixel filtering with L8 Radiance filter
% DRad_L8radfilt = D_TOArad_final./L8_Rad_image_mask;
% figure, imagesc(DRad_L8radfilt); colorbar
% figure, histogram(DRad_L8radfilt)

% Dove Rad pixel filtering with L8 SNR filter
DRad_L8SNRfilt = D_TOArad_final./L8_SNR_image_mask;
% figure, imagesc(DRad_L8SNRfilt); colorbar
% figure, histogram(DRad_L8SNRfilt)

% Dove SD pixel filtering with L8 SNR filter
D_SD_L8SNRfilt = D_SD_image./L8_SNR_image_mask;
% figure, imagesc(D_SD_L8SNRfilt); colorbar
% figure, histogram(D_SD_L8SNRfilt)

% Dove SNR
Dove_SNR = DRad_L8SNRfilt./D_SD_L8SNRfilt;
% figure, imagesc(Dove_SNR); colorbar
% figure, histogram(Dove_SNR)

DoveL8SNR.L8SNR{scene, :} = L8_SNR_image_filtered;
DoveL8SNR.DoveSNR{scene, :} = Dove_SNR;

Dove_max_SNR(scene) = max(max(Dove_SNR));
Dove_min_SNR(scene) = min(min(Dove_SNR));

L8_max_SNR
L8_min_SNR

Dove_max_SNR
Dove_min_SNR
 end

%% Saving the data
AllSNR.PR038037.L8SNR = DoveL8SNR.L8SNR;
AllSNR.PR038037.DoveSNR = DoveL8SNR.DoveSNR;

AllSNR.PR038038.L8SNR = DoveL8SNR.L8SNR;
AllSNR.PR038038.DoveSNR = DoveL8SNR.DoveSNR;

%%% The following 3 Sections are for checking the Dove and L8 SNR
%%% Path/Row - 038/037- Scene- 16, 22, 29
%% Scene 16
DSNR_S16 = AllSNR.PR037038.DoveSNR{1, 1};
% figure, imagesc(DSNR_S16); colorbar
% figure, histogram(DSNR_S16, 100)
DSNR_S16_f = DSNR_S16(find(~isnan(DSNR_S16(:))));

L8SNR_S16 = AllSNR.PR037038.L8SNR{1, 1};
% figure, imagesc(L8SNR_S16); colorbar
% figure, histogram(L8SNR_S16, 100)
L8SNR_S16_f = L8SNR_S16(find(~isnan(L8SNR_S16(:))));
%figure, plot(L8SNR_S16_f, DSNR_S16_f, '*')

L8SNR_S16_filtered = L8SNR_S16_f(L8SNR_S16_f>360 & L8SNR_S16_f< 370);
DSNR_S16_filtered = DSNR_S16_f(L8SNR_S16_f>360 & L8SNR_S16_f< 370);

%% Scene 22
DSNR_S22 = AllSNR.PR037038.DoveSNR{2, 1};
% figure, imagesc(DSNR_Scene22); colorbar
% figure, histogram(DSNR_Scene22, 100)
% [Y X] = find(~isnan(DSNR_Scene22));
DSNR_S22_f = DSNR_S22(find(~isnan(DSNR_S22(:))));

L8SNR_S22 = AllSNR.PR037038.L8SNR{2, 1};
% figure, imagesc(L8SNR_S22); colorbar
% figure, histogram(L8SNR_S22, 100)
L8SNR_S22_f = L8SNR_S22(find(~isnan(L8SNR_S22(:))));
% figure, plot(L8SNR_S22_f, DSNR_S22_f, '*')

L8SNR_S22_filtered = L8SNR_S22_f(L8SNR_S22_f>360 & L8SNR_S22_f< 370)
L8SNR_S22_filtered = DSNR_S22_f(L8SNR_S22_f>360 & L8SNR_S22_f< 370)

%% Scene 29
DSNR_S29 = AllSNR.PR037038.DoveSNR{3, 1};
% figure, imagesc(DSNR_Scene29); colorbar
% figure, histogram(DSNR_Scene29, 100)
DSNR_S29_f = DSNR_S29(find(~isnan(DSNR_S29(:))));
 
L8SNR_S29 = AllSNR.PR037038.L8SNR{3, 1};
% figure, imagesc(L8SNR_S29); colorbar
% figure, histogram(L8SNR_S29, 100)
L8SNR_S29_f = L8SNR_S29(find(~isnan(L8SNR_S29(:))));
% figure, plot(L8SNR_S29_f, DSNR_S29_f, '*')

L8SNR_S29_filtered = L8SNR_S29_f(L8SNR_S29_f>360 & L8SNR_S29_f< 370)
DSNR_S29_filtered = DSNR_S29_f(L8SNR_S29_f>360 & L8SNR_S29_f< 370)

%% The following 4 sections are for 
%%% Path/Row- 038/038- Scene 36, 42, 49, 55
%% Scene 36
DSNR_S36 = AllSNR.PR038038.DoveSNR{1, 1};
% figure, imagesc(DSNR_S36); colorbar
% figure, histogram(DSNR_S36, 100)
DSNR_S36_f = DSNR_S36(find(~isnan(DSNR_S36(:))));

L8SNR_S36 = AllSNR.PR038038.L8SNR{1, 1};
% figure, imagesc(L8SNR_S36); colorbar
% figure, histogram(L8SNR_S36, 100)
L8SNR_S36_f = L8SNR_S36(find(~isnan(L8SNR_S36(:))));
%figure, plot(L8SNR_S36_f, DSNR_S36_f, '*')

L8SNR_S36_filtered = L8SNR_S36_f(L8SNR_S36_f>360 & L8SNR_S36_f< 370)
DSNR_S36_filtered = DSNR_S36_f(L8SNR_S36_f>360 & L8SNR_S36_f< 370)

%% Scene 42
DSNR_S49 = AllSNR.PR038038.DoveSNR{2, 1};
% figure, imagesc(DSNR_S42); colorbar
% figure, histogram(DSNR_S42, 100)
DSNR_S42_f = DSNR_S49(find(~isnan(DSNR_S49(:))));

L8SNR_S42 = AllSNR.PR038038.L8SNR{2, 1};
% figure, imagesc(L8SNR_S42); colorbar
% figure, histogram(L8SNR_S42, 100)
L8SNR_S42_f = L8SNR_S42(find(~isnan(L8SNR_S42(:))));
%figure, plot(L8SNR_S42_f, DSNR_S42_f, '*')

L8SNR_S42_filtered = L8SNR_S42_f(L8SNR_S42_f>360 & L8SNR_S42_f< 370)
DSNR_S42_filtered = DSNR_S42_f(L8SNR_S42_f>360 & L8SNR_S42_f< 370)
% figure, histogram(L8SNR_S42_filtered, 100)
% figure, histogram(DSNR_S42_filtered, 100)

%% Scene 49
DSNR_S49 = AllSNR.PR038038.DoveSNR{3, 1};
% figure, imagesc(DSNR_S49); colorbar
% figure, histogram(DSNR_S49, 100)
DSNR_S49_f = DSNR_S49(find(~isnan(DSNR_S49(:))));

L8SNR_S49 = AllSNR.PR038038.L8SNR{2, 1};
% figure, imagesc(L8SNR_S49); colorbar
% figure, histogram(L8SNR_S49, 100)
L8SNR_S49_f = L8SNR_S49(find(~isnan(L8SNR_S49(:))));
%figure, plot(L8SNR_S42_f, DSNR_S42_f, '*')

L8SNR_S49_filtered = L8SNR_S49_f(L8SNR_S49_f>360 & L8SNR_S49_f< 370)
DSNR_S49_filtered = DSNR_S49_f(L8SNR_S49_f>360 & L8SNR_S49_f< 370)

%% Scene 55
DSNR_S49 = AllSNR.PR038038.DoveSNR{3, 1};
% figure, imagesc(DSNR_S49); colorbar
% figure, histogram(DSNR_S49, 100)
DSNR_S49_f = DSNR_S49(find(~isnan(DSNR_S49(:))));

L8SNR_S49 = AllSNR.PR038038.L8SNR{2, 1};
% figure, imagesc(L8SNR_S49); colorbar
% figure, histogram(L8SNR_S49, 100)
L8SNR_S49_f = L8SNR_S49(find(~isnan(L8SNR_S49(:))));
%figure, plot(L8SNR_S42_f, DSNR_S42_f, '*')

L8SNR_S49_filtered = L8SNR_S49_f(L8SNR_S49_f>360 & L8SNR_S49_f< 370)
DSNR_S49_filtered = DSNR_S49_f(L8SNR_S49_f>360 & L8SNR_S49_f< 370)