clear all

base='Z:\ImageDrive\PlanetLabs\Processed\1047\P039\R037\'; 
dates = dir(base); dates([1 2])=[]; 

% Dove File Path
DoveDataPath = fullfile(base, dates.name, '*_AnalyticMS.tif');
DoveDirectory = dir(DoveDataPath);
% DoveDataPathXML = fullfile(base, dates.name, '*_metadata.xml');
DoveDataPathXML = fullfile(base, dates.name, '*_metadata.xml');
DoveDirectoryXML = dir(DoveDataPathXML);
[NoOf_Scene, y1] = size(DoveDirectory);

for scene = 3 
    for D_band = 1 %1:4 band = 2
        %%% Dove
        % Dove Image file- scene by scene
        Dove_image_file = fullfile(base, dates.name, DoveDirectory(scene).name);

        % Reading Dove Image
        Dove_image_all_band =imread(Dove_image_file);
        DN_Dove = Dove_image_all_band(:,:, D_band);
        DN_Dove = double(DN_Dove);
        DN_Dove(DN_Dove == 0) = NaN;
        band_sf = 0.01;
        Dove_TOArad = DN_Dove*band_sf;
        
        %%%% Looking at the image
        % figure, imagesc(Dove_TOArad); colorbar
        % figure, histogram(Dove_TOArad)
        

    end
    
  R_dove = [0 -3; 3 0; 682731 3612660];
  map2pix(R_dove, 528570, 3601680);

end 
