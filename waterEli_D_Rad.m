function [Dove_TOArad_w_filter] = waterEli_D_Rad(base, D_band, scene_no, threshold)

%% Water Elimination Filter for DOve
dates = dir(base); dates([1 2])=[]; 
D_band = D_band; scene = scene_no;
threshold = threshold;

% Dove Image File Path
DoveDataPath = fullfile(base, dates.name, '*_AnalyticMS.tif');
DoveDirectory = dir(DoveDataPath);

% Dove XML File Path
DoveDataPathXML = fullfile(base, dates.name, '*_metadata.xml');
DoveDirectoryXML = dir(DoveDataPathXML); [NoOf_Scene, y] = size(DoveDirectory);

% Extracting all the values from Metadata file
MData_file = fullfile(base, dates.name, DoveDirectoryXML(scene).name);
[MData_values]= xml2struct_new_v(MData_file);
[MData_values_old]= xml2struct(MData_file);

%%% Dove
% Dove Image file- scene by scene
Dove_image_file = fullfile(base, dates.name, DoveDirectory(scene).name);

% Reading Dove Image
Dove_image_all_band = imread(Dove_image_file);
DN_Dove = Dove_image_all_band(:,:, D_band);
DN_Dove = double(DN_Dove);
DN_Dove(DN_Dove == 0) = NaN;
band_sf = 0.01;
Dove_TOArad = DN_Dove*band_sf;

%%%% Looking at the image
% figure, imagesc(Dove_TOArad); colorbar
% figure, histogram(Dove_TOArad)

%%%% Sun Elevation Angles and Azimuth Angles 
Dove_image_info.D_SEle_Angle = str2num(MData_values.ps_colon_EarthObservation.gml_colon_using...
            .eop_colon_EarthObservationEquipment.eop_colon_acquisitionParameters.ps_colon_Acquisition...
            .opt_colon_illuminationElevationAngle.Text);  
Dove_image_info.D_SAzi_Angle=str2num(MData_values.ps_colon_EarthObservation.gml_colon_using...
           .eop_colon_EarthObservationEquipment.eop_colon_acquisitionParameters.ps_colon_Acquisition...
           .opt_colon_illuminationAzimuthAngle.Text);  
Dove_TOArad = Dove_TOArad./sind(Dove_image_info.D_SEle_Angle); % Cosine Correction
%%%% Looking at the image
% figure, imagesc(Dove_TOArad); colorbar
% figure, histogram(Dove_TOArad)

%%%% Filtering out the low radiance data 
Dove_TOArad_temp = Dove_TOArad;
Dove_TOArad_temp(Dove_TOArad_temp < threshold) = nan;
%%%% Looking at the image
% figure, imagesc(Dove_TOArad_temp); colorbar
% figure, histogram(Dove_TOArad_temp)

%%%% Creating water mask
Dove_TOArad_w_filter = Dove_TOArad_temp;
Dove_TOArad_w_filter(~isnan(Dove_TOArad_w_filter)) = 1;
% figure, imagesc(Dove_TOArad_w_filter); colorbar
end
