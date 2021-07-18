clear all
%%% Directory 
% base='Z:\SpecialNeeds\Sakib\DoveData\P039\R037';
base='Z:\ImageDrive\PlanetLabs\Processed\DoveData\dove_r_2171395_1058\scenes';
dates = dir(base); dates([1 2])=[];

% Dove File Path
DoveDataPath = fullfile(base, '*_AnalyticMS.tif');
DoveDirectory = dir(DoveDataPath);
DoveDataPathXML = fullfile(base, '*_metadata.xml');
DoveDirectoryXML = dir(DoveDataPathXML);
[NoOf_Scene, y] = size(DoveDirectory);

%%
scene = 3;
  for window_D = 3:2:61
   for D_band=1 %:4
    % Dove Image file- scene by scene
    Dove_image_file = fullfile(base, DoveDirectory(scene).name);
    
    % Reading Dove Image
    Dove_image_all_band =imread(Dove_image_file);
    DN_Dove = Dove_image_all_band(:,:, D_band);
    DN_Dove = double(DN_Dove);
    DN_Dove(DN_Dove == 0) = NaN;
    band_sf = 0.01;
    Dove_TOArad = DN_Dove*band_sf;
 
    % Mean and Standard Deviation Image
    Dove_Lmean_image = conv2(Dove_TOArad, ones(window_D)./(window_D*window_D),'same'); 
    k = (window_D-1)/2;  % window = 2*k+1; k is defined in this way in movingstd2 
    D_LSD_image = movingstd2(Dove_TOArad, k);

    % SNR calculation
    SNR_esti_D = Dove_Lmean_image./D_LSD_image;
    
    % Checking the SNR image and SNR distribution
    % figure; histogram(D_LSD_image)
    
    %%% Bin Width Selection
    LSD_min = min(min(D_LSD_image));
    LSD_Avg = nanmean(nanmean(D_LSD_image));
    % figure, histogram(D_LSD_image, 150, 'BinLimits', [LSD_min 1.2*LSD_Avg])
    
    % Storing the Noise at Maximum box count
    [N, edges] = histcounts(D_LSD_image,  150, 'BinLimits',[LSD_min 1.2*LSD_Avg]);
    maximum= max(max(N));
    [x, y]= find (N== maximum);
    Noise(window_D)= edges(y);
    % max_count(window) = max(N);
    window_size(window_D)= window_D;
%         
%     Noise(Noise==0)= NaN;
%     N(N<=200)= NaN;
%     minimum= min(min(N));
%     [i, j]=find(N==minimum);
%     Noise_highest(window_D) = edges(max(j));
%     Noise_highest(Noise_highest==0)= NaN;
   end  
  
  end 
% figure, plot(D_SD_image, Dove_mean_image, 'c.')
% figure, plot(window_size, Noise, '*')
