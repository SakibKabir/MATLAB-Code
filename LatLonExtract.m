clear all
% base location
base='Z:\ImageDrive\PlanetLabs\Processed\1047\P039\R037\';
dates = dir(base);
dates([1 2])=[];
% Xml file 
MData_file=strcat(base, dates.name, filesep, '20181223_172655_38_1047_3B_AnalyticMS_metadata.xml');
[MData_values] = xml2struct_new_v(MData_file);

% Latitude and Longitude
lat_TL = str2num(MData_values.ps_colon_EarthObservation.gml_colon_target...
    .ps_colon_Footprint.ps_colon_geographicLocation.ps_colon_topLeft.ps_colon_latitude.Text);   
lon_TL = str2num(MData_values.ps_colon_EarthObservation.gml_colon_target...
    .ps_colon_Footprint.ps_colon_geographicLocation.ps_colon_topLeft.ps_colon_longitude.Text);
