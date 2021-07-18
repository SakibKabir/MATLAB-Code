function [Dove_TOAref_w_filter] = waterEli_D_Ref(base_D, D_band, scene_no, percentage)
    %% Water Elimination Filter for Dove TOA Reflectance
    dates = dir(base_D); dates([1 2])=[]; 
    prc = percentage;
    
    % Dove Image File Path
    DoveDataPath = fullfile(base_D, dates.name, '*_AnalyticMS.tif');
    DoveDirectory = dir(DoveDataPath)

    % Dove XML File Path
    DoveDataPathXML = fullfile(base_D, dates.name, '*_metadata.xml')
    DoveDirectoryXML = dir(DoveDataPathXML);
    
    % Extracting all the values from Metadata file
    MData_file = fullfile(base_D, dates.name, DoveDirectoryXML(scene_no).name)
    [MData_values]= xml2struct_new_v(MData_file);
    [MData_values_old]= xml2struct(MData_file);

    % Dove Image file- scene by scene
    Dove_image_file = fullfile(base_D, dates.name, DoveDirectory(scene_no).name);

    % Reading Dove Image
    Dove_image_all_band = imread(Dove_image_file);
    DN_Dove = Dove_image_all_band(:,:, D_band);
    DN_Dove = double(DN_Dove);
    DN_Dove(DN_Dove == 0) = NaN;
    
    %%% Reflectance conversion coefficient
    if D_band == 1
        id = 12; 
    elseif D_band == 2
        id = 14;
    elseif D_band == 3
        id = 16;
    elseif D_band == 4
        id = 18;
    end
    
    %%%% Dove TOA Reflectance Conversion
    Ref_coeff = str2double(MData_values_old.Children(10).Children(2).Children(id).Children(10).Children.Data);
    Dove_TOAref = DN_Dove*Ref_coeff;
    
    % Landsat 8 water filter from NIR band
    Dove_TOAref_prc = Dove_TOAref;
    Dove_TOAref_prc = Dove_TOAref_prc(~isnan(Dove_TOAref_prc));

    Dove_TOAref_prc_sorted = sort(Dove_TOAref_prc);
    ind_prc = round(prc*length(Dove_TOAref_prc_sorted)); % getting the index at prc
    Data_prc = Dove_TOAref_prc_sorted(find(Dove_TOAref_prc_sorted) == ind_prc); % getting data at prc

    % Filtering out the data
    Dove_TOAref_filtered = Dove_TOAref;
    Dove_TOAref_filtered(Dove_TOAref_filtered <= Data_prc) = NaN; 
    % figure, imagesc(L8_TOAref_filtered); colorbar
    
    %%%% Creating water mask
    Dove_TOAref_w_filter = Dove_TOAref_filtered;
    Dove_TOAref_w_filter(~isnan(Dove_TOAref_w_filter)) = 1;
    % figure, imagesc(Dove_TOAref_w_filter); colorbar
end
