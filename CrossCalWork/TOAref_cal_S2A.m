function [S2A_TOAref, R_S2A] = TOAref_cal_S2A(base_S2A, date_S2A, band)
    %%% Sentinel 2A: TOA Reflectance and Reference Matrix
    % Here 8A band is band 13
    MSI2A_bands =  {'01' '02' '03' '04' '05' '06' '07' '08' '09' '10' '11' '12' '8A'};
    % Resolution =  [60   10   10   10   20   20   20   10   60   60   20   20   20];
 
    %%% 2 XML File Parsing 
    MSIpath = dir(fullfile(base_S2A, date_S2A, 'MSIL1C', '*MSIL1C.xml'));
    MTLpath = dir(fullfile(base_S2A, date_S2A, 'MSIL1C', '*TL.xml'));
    MSIfile=fullfile(base_S2A, date_S2A, 'MSIL1C', MSIpath.name);
    MTDfile=fullfile(base_S2A, date_S2A, 'MSIL1C', MTLpath.name);
    % Storing the values from xml file
    [MSI_values]= xml2struct_new_v(MSIfile);
    [MTD_TL] = xml2struct_new_v(MTDfile);

    % Image file location
    Band_No = dir(fullfile(base_S2A, date_S2A, 'MSIL1C', strcat('*B', MSI2A_bands{band}, '.jp2')));
    S2_image_file = fullfile(base_S2A, date_S2A, 'MSIL1C', Band_No.name);

    % Reading S2A image
    DN_S2A = imread(S2_image_file);     
    DN_S2A = double(DN_S2A);
    DN_S2A(DN_S2A == 0)= nan;
    S2A_info = imfinfo(S2_image_file);
    % S2A TOA Reflectance
    Q_value = 10000;
    S2A_TOAref = DN_S2A./Q_value;
    % figure, imagesc(S2A_TOAref); colorbar;

    %%
    %Extracting the necsseary values from xml file
    %BandID selection depending on band number
    if band == 2 || band ==3 ||band == 4 || band == 8
        BandID = 1; %% 10m

    elseif band == 1 ||band == 9 ||band == 10
        BandID = 3; %% 60m

    elseif  band == 5 || band == 6 ||band == 7 ||band == 11 ||band == 12 ||band == 13
        BandID = 2; %% 20m
    end

    %%
    s2_image_info.CornerCoords.X = str2num(MTD_TL.n1_colon_Level_dash_1C_Tile_ID....
            .n1_colon_Geometric_Info.Tile_Geocoding.Geoposition{1, BandID}.ULX.Text);

    s2_image_info.CornerCoords.Y = str2num(MTD_TL.n1_colon_Level_dash_1C_Tile_ID....
            .n1_colon_Geometric_Info.Tile_Geocoding.Geoposition{1, BandID}.ULY .Text);

    s2_image_info.Height = str2num(MTD_TL.n1_colon_Level_dash_1C_Tile_ID....
            .n1_colon_Geometric_Info.Tile_Geocoding.Size{1, BandID}.NROWS.Text);

    s2_image_info.Width = str2num(MTD_TL.n1_colon_Level_dash_1C_Tile_ID....
            .n1_colon_Geometric_Info.Tile_Geocoding.Size{1, BandID}.NCOLS.Text);

    s2_image_info.CurrentRes = str2num(MTD_TL.n1_colon_Level_dash_1C_Tile_ID.n1_colon_Geometric_Info....
            .Tile_Geocoding.Geoposition{1, BandID}.Attributes.resolution);

    %Reference Matrix 
    R_S2A = [0 -s2_image_info.CurrentRes; s2_image_info.CurrentRes 0;...
             s2_image_info.CornerCoords.X  s2_image_info.CornerCoords.Y];
end