clear all
%%% Directory 
base='Z:\ImageDrive\OLI.TIRS\L8\P023\R036';
dates = dir(base); dates([1 2])=[];
bands=  {'1' '2' '3' '4' '5' '6' '7' '8'}; 

%MTL File Parsing 
MTL=dir(fullfile(base, dates.name,'*MTL.txt'));
[MTL_List_L8, vaule]= MTL_parser_L8(fullfile(MTL.folder, MTL.name));

for L8_band = 2 %:4
    L8_Folder_info=dir(fullfile(base, dates.name, strcat('*B', bands{L8_band},'.tif')));  
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
end

    %%% OLI Image 
    % figure, imagesc(L8_TOArad ); colorbar
    L8_TOArad_filtered = L8_TOArad; % for filtering TOA Rad with typical
    %%% OLI Radiance Image filtering
    L8_TOArad_filtered(L8_TOArad_filtered < 39.5 | L8_TOArad_filtered > 40) = NaN;
    % figure, imagesc(L8_TOArad_filtered); colorbar

    %%% Dove Standard Deviation image
    window = 5;
    k= (window-1)/2; 
    L8_SD_image = movingstd2(L8_TOArad, k);
    % figure, imagesc(L8_SD_image); colorbar

    %%% SNR 
    L8_SNR_filtered = L8_TOArad_filtered./L8_SD_image;
    % figure, imagesc(L8_SNR_filtered); colorbar
    % figure, histogram(L8_SNR_filtered)
    % L8SNR.allSNR{scene, :} = L8_SNR_filtered;
    
    %%% Meand and SD calculation
    temp_L8_SNR_filtered = L8_SNR_filtered;
    temp_L8_SNR_filtered = temp_L8_SNR_filtered(~isnan(temp_L8_SNR_filtered));
    Mean = round(mean(temp_L8_SNR_filtered), 1);
    Std_L8 = std(temp_L8_SNR_filtered);
    
    SDs = round((Mean + (1:6)*Std_L8),1);
    %figure(scene),
    h = histogram(L8_SNR_filtered); % Histogram
    N = max(h.Values); % Maximum bin count

    hold on
    plot([Mean Mean],[0 N],'r','LineWidth',3) % Mean
    X = repmat(Mean+(1:6)*Std_L8, 2, 1);
    Y = repmat([0;N],1,6);
    plot(X,Y,'Color',[0 0 0],'LineWidth',3) % Standard deviations
    legend('Data',['Mean (' num2str(Mean) ')'], ['SDs[1-6] (' num2str(SDs) ')' ])
    xlabel('Signal to Noise Ratio')
    ylabel('Number of Observation')
    ax  = gca;
    ax.FontSize = 36;
    hold off