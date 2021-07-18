clear all
%%% Directory 
base='Z:\ImageDrive\OLI.TIRS\L8\P023\R033\';
dates = dir(base); dates([1 2])=[];
bands=  {'1' '2' '3' '4' '5' '6' '7' '8'}; 
%MTL File Parsing 
MTL = dir(fullfile(base, dates.name, 'LC1','*MTL.txt'));
[MTL_List_L8, vaule]= MTL_parser_L8(fullfile(MTL.folder, MTL.name));

%%
 for window_L8 = 3:2:31 %:51%61
   for L8_band = 2 %:4
    L8_Folder_info=dir(fullfile(base, dates.name,'LC1', strcat('*B', bands{L8_band},'.tif')));  
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
    
    % Mean and Standard Deviation Image
    L8_mean_image = conv2(L8_TOArad, ones(window_L8)./(window_L8*window_L8),'same'); 
    k = (window_L8-1)/2;  % window = 2*k+1; k is defined in this way in movingstd2 
    L8_SD_image = movingstd2(L8_TOArad, k);
    %Mean_SD(window) = nanmean(nanmean(Dove_SD_image));
        
    % SNR calculation
    SNR_esti_L8 = L8_mean_image./L8_SD_image;
    
    % Checking the SNR image and SNR distribution
    % figure; histogram(L8_SD_image)
    
    [N, edges] = histcounts(L8_SD_image);
    maximum= max(max(N));
    [x, y]= find (N== maximum);
    Noise(window_L8)= edges(y);
    % max_count(window) = max(N);
    window_size(window_L8)= window_L8;
        
    Noise(Noise==0)= NaN;
    N(N<=200)= NaN;
    minimum= min(min(N));
    [i, j]=find(N==minimum);
    Noise_highest(window_L8) = edges(max(j));
    Noise_highest(Noise_highest==0)= NaN;

  end  
    
 end 
% figure, plot(L8_SD_image, L8_mean_image, 'c.')
% figure, plot(window_size, Noise, '*')