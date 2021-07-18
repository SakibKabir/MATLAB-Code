%% Landsat 8 TOA Reflectance Comparision with vegetation hyperspectral (SDSU Cal Site data) 
% Storing all the TOA reflectance band by band
function [L8_allMean, L8_allSD, D_allMean, D_allSD] = OLI_DoveTOACCdata()
%clear all
%load('CCResults_unstructuredfull2.mat')
load('CCTOAref.mat')
%%
for location = 1:19
    noofScene = size(CCTOAref.L8_mean{location, 1},1);
    for scene_no = 1:noofScene
        L8_refbands = CCTOAref.L8_mean{location, 1}{scene_no, 1};
        L8_sdbands = CCTOAref.L8_sd{location, 1}{scene_no, 1};
        
        D_refbands = CCTOAref.D_mean{location, 1}{scene_no, 1}; 
        D_sdbands = CCTOAref.D_sd{location, 1}{scene_no, 1}; 
                  
        % Filtering all the fill values using Dove as mask 
        L8_refbands = L8_refbands.*(~isnan(D_refbands));
        L8_refbands(L8_refbands ==0) = nan;
        L8_refbands(any(isnan(L8_refbands), 2), :) = [];

        L8_sdbands = L8_sdbands.*(~isnan(D_refbands));
        L8_sdbands(L8_sdbands ==0) = nan;
        L8_sdbands(any(isnan(L8_sdbands), 2), :) = [];
        
        % Removing all the nan's from Dove Matrix
        D_refbands(any(isnan(D_refbands), 2), :) = [];
        D_sdbands(any(isnan(D_sdbands), 2), :) = [];
        
        % Storing data scene by scene
        L8_SbyS(scene_no,:) = {L8_refbands}; clear L8_refbands
        L8_sdSbyS(scene_no,:) = {L8_sdbands}; clear L8_sdbands
        D_SbyS(scene_no,:) = {D_refbands}; clear D_refbands
        D_sdSbyS(scene_no,:) = {D_sdbands}; clear D_sdbands
        
        % Merging all Scenes data in one variable
        L8_allSmean = cell2mat(L8_SbyS); 
        L8_allSsd = cell2mat(L8_sdSbyS); 
        D_allSmean = cell2mat(D_SbyS);
        D_allSsd = cell2mat(D_sdSbyS);
        
    end; clear L8_SbyS L8_sdSbyS D_SbyS D_sdSbyS
    
        % Storing all data location by location
        L8_LbyL(location, :) = {L8_allSmean}; %clear L8_allSmean
        L8sd_LbyL(location, :) = {L8_allSsd};% clear L8_allSsd
        D_LbyL(location, :) = {D_allSmean}; %clear D_allSmean
        Dsd_LbyL(location, :) = {D_allSsd};% clear D_allSsd
end

L8_LbyL = L8_LbyL(~cellfun('isempty',L8_LbyL));
L8sd_LbyL = L8sd_LbyL(~cellfun('isempty',L8sd_LbyL));
D_LbyL = D_LbyL(~cellfun('isempty',D_LbyL));
Dsd_LbyL = Dsd_LbyL(~cellfun('isempty',Dsd_LbyL));

No_L = size(L8_LbyL, 1);

%% Structuring all the data for plotting
for location = 1:No_L
    L8_meanTOAref = L8_LbyL{location, 1}; L8_sdTOAref = L8sd_LbyL{location, 1};  
    D_meanTOAref =  D_LbyL{location, 1}; D_sdTOAref =  Dsd_LbyL{location, 1};

    L8_Blue(location,:) ={L8_meanTOAref(:,1)}; L8_sdBlue(location,:) ={L8_sdTOAref(:,1)};    
    D_Blue(location,:) ={D_meanTOAref(:,1)}; D_sdBlue(location,:) ={D_sdTOAref(:,1)};
    
    L8_Green(location,:) ={L8_meanTOAref(:,2)}; L8_sdGreen(location,:) ={L8_sdTOAref(:,2)};
    D_Green(location,:) ={D_meanTOAref(:,2)}; D_sdGreen(location,:) ={D_sdTOAref(:,2)};
    
    L8_Red(location,:) ={L8_meanTOAref(:,3)}; L8_sdRed(location,:) ={L8_sdTOAref(:,3)};  
    D_Red(location,:) ={D_meanTOAref(:,3)}; D_sdRed(location,:) ={D_sdTOAref(:,3)};
    
    L8_NIR(location,:) ={L8_meanTOAref(:,4)}; L8_sdNIR(location,:) ={L8_sdTOAref(:,4)};  
    D_NIR(location,:) ={D_meanTOAref(:,4)}; D_sdNIR(location,:) ={D_sdTOAref(:,4)};
end; %clear L8_meanTOAref L8_sdTOAref D_meanTOAref D_sdTOAref

%% Band by band All data 
L8_Blue_all = cell2mat(L8_Blue); L8_sdBlue_all = cell2mat(L8_sdBlue); D_Blue_all = cell2mat(D_Blue); D_sdBlue_all = cell2mat(D_sdBlue);
L8_Green_all = cell2mat(L8_Green); L8_sdGreen_all = cell2mat(L8_sdGreen); D_Green_all = cell2mat(D_Green); D_sdGreen_all = cell2mat(D_sdGreen);
L8_Red_all = cell2mat(L8_Red); L8_sdRed_all = cell2mat(L8_sdRed); D_Red_all = cell2mat(D_Red); D_sdRed_all = cell2mat(D_sdRed);
L8_NIR_all = cell2mat(L8_NIR); L8_sdNIR_all = cell2mat(L8_sdNIR); D_NIR_all = cell2mat(D_NIR); D_sdNIR_all = cell2mat(D_sdNIR);


L8_allMean = [cell2mat(L8_Blue) cell2mat(L8_Green) cell2mat(L8_Red) cell2mat(L8_NIR)]; 
L8_allSD = [cell2mat(L8_sdBlue) cell2mat(L8_sdGreen) cell2mat(L8_sdRed) cell2mat(L8_sdNIR)]; 
D_allMean = [cell2mat(D_Blue) cell2mat(D_Green) cell2mat(D_Red) cell2mat(D_NIR)]; 
D_allSD = [cell2mat(D_sdBlue) cell2mat(D_sdGreen) cell2mat(D_sdRed) cell2mat(D_sdNIR)]; 
%figure, plot(L8_allMean(:,1),D_allMean(:,1), '.' )
end
