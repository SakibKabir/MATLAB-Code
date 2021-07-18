function [Dove_TOAref, R_Dove, R_doveEx] = TOAref_cal_D(base_D, scene_no, D_band, equator)
    %% Dove TOA Reflectance Calculation Function
    dates = dir(base_D); 
    dates([1 2])=[];
    
    % Dove Image File Path
    DoveDataPath = fullfile(base_D, dates.name, '*_AnalyticMS_reprojected.tif');
    DoveDirectory = dir(DoveDataPath);

    % Dove XML File Path
    DoveDataPathXML = fullfile(base_D, dates.name, '*_metadata.xml');
    DoveDirectoryXML = dir(DoveDataPathXML);
    [NoOf_Scene, y] = size(DoveDirectory);

    % Extracting all the values from Metadata file
    MData_file = fullfile(base_D, dates.name, DoveDirectoryXML(scene_no).name);
    [MData_values]= xml2struct_new_v(MData_file);
    [MData_values_old]= xml2struct(MData_file);
    
    % Dove Image file- scene by scene
    Dove_image_file = fullfile(base_D, dates.name, DoveDirectory(scene_no).name);

    % Geting the Dove Reference matrix
    Image_info = imfinfo(Dove_image_file);
    ModelTiepoinTag = Image_info.ModelTiepointTag;
    X_TL = ModelTiepoinTag(4);
    Y_TL = ModelTiepoinTag(5);
    R_Dove = [0 -3; 3 0; X_TL Y_TL-equator];
    
    % Reading Dove Image
    Dove_image_all_band = imread(Dove_image_file);
    [Dove, R_doveEx] = geotiffread(Dove_image_file);
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
    
    %%%% Looking at the image and its distribution
    % figure, imagesc(Dove_TOAref); colorbar
    % figure, histogram(Dove_TOAref)
end