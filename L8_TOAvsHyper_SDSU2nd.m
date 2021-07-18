clear all
%% Loading the Hyperspectral profiles and OLI Reflectance
[Hyper_WLenNM, Hyper_Spectrum, Hy_WL_Sand, Hyper_Spectrum_all] = HyperData();
L8_allMean = OLI_DoveTOACCdata(); 

%% Calculating SBAF with Appropiate Profile
for ref = 1%:size(L8_allMean) % all the reflectances inside the ROI
       OLI_ref = L8_allMean(ref,:); % OLI Reflectance ROI by ROI
 for NoPro = 1:size(Hyper_Spectrum_all, 2)
        %% For Vegitation
        Hy_Spectrum_Vegi = Hyper_Spectrum_all(1:size(Hyper_Spectrum,1), NoPro);
        %%% Banding Reference Sensor: OLI- Satnum = 24
        [Labels_L8_Vegi, HyBanded_L8_Vegi] = bander(Hyper_WLenNM, Hy_Spectrum_Vegi, 24); clear Labels_L8_Vegi
        HyB_L8 = HyBanded_L8_Vegi(2:5); clear Labels_L8_Vegi HyBanded_L8_Vegi
        M_vegi = mean(OLI_ref./HyB_L8);
        
        if M_vegi > 1 % Means OLI Ref is higher; so, multiply Hyper Ref to scale it up
            OLI_ref_m = OLI_ref;
            HyB_L8_m = HyB_L8.*M_vegi;
        elseif M_vegi < 1 % Means OLI Ref is lower; so, multiply OLI Ref to scale it up
            OLI_ref_m = OLI_ref.*M_vegi;
            HyB_L8_m = HyB_L8;
        end; clear HyB_L8
        
        if isnan(M_vegi)
            SSE(1, NoPro) = nan;
        else
            SSE(1, NoPro) = sse(OLI_ref_m, HyB_L8_m);
        end; clear M_vegi Hy_Spectrum_Vegi HyBanded_L8_Vegi OLI_ref_m HyB_L8_m
        
       if NoPro < 24
        %%% For Sand
        Hy_Spectrum_Sand = Hyper_Spectrum_all(size(Hyper_Spectrum, 1)+1:end, NoPro);
        %%% Banding Reference Sensor: OLI- Satnum = 24
        [Labels_L8_Sand, HyBanded_L8_sand] = bander(Hy_WL_Sand, Hy_Spectrum_Sand, 24); clear Hy_Spectrum_Sand Labels_L8_Sand
        HyB_L8s = HyBanded_L8_sand(2:5); clear HyBanded_L8_sand 
        M_sand = mean(OLI_ref./HyB_L8s);
        
        if M_sand > 1
            OLI_ref_m = OLI_ref;
            HyB_L8_m = HyB_L8s.*M_sand;
        elseif M_sand < 1
            OLI_ref_m = OLI_ref.*M_sand;
            HyB_L8_m = HyB_L8s;
        end; clear M_sand HyB_L8s
        
        % I Know I don't have any nan for sand
        SSE(2, NoPro) = sse(OLI_ref_m, HyB_L8_m);  
       else
        %clear Labels_L8_Sand HyBanded_L8_sand OLI_ref_m HyB_L8_m
       end 
 end; clear NoPro       
    
SSE(SSE==0) = nan; [r, c]= find(SSE == min(min(SSE))); %clear SSE
%%
if r(1) == 1 % Vegi
    Hyper = Hyper_Spectrum_all(1:size(Hyper_Spectrum, 1), c(1));
    Hyper_WL = Hyper_WLenNM; 

elseif r(1) == 2 % Sand
    Hyper = Hyper_Spectrum_all(size(Hyper_Spectrum, 1)+1:end, c(1));
    Hyper_WL = Hy_WL_Sand;
end; clear r c 

    %%% Banding Reference Sensor OLI- Satnum = 24
    [labels_L8, HyBanded_L8] = bander(Hyper_WL, Hyper, 24); clear labels_L8 
    HyBanded_OLIRef = HyBanded_L8(:, 2:5); clear HyBanded_L8
    %%% Banding Target Sensor Dove R(1047)- Satnum = 84
    [labels_1047, HyBanded_1047] = bander(Hyper_WL, Hyper, 84);  clear labels_1047 Hyper_WL Hyper 
    SBAF_L8_D1047_SCalSite(ref, :) = HyBanded_OLIRef./HyBanded_1047; clear HyBanded_1047
    
    %%% Store the SSE of OLI and Hyper Banded Reflectance
     M = mean(OLI_ref./HyBanded_OLIRef);
        if M > 1 % Means OLI Ref is higer; so, multiply Hyper Ref to scale it up
            OLI_ref_m = OLI_ref;
            HyB_L8_m = HyBanded_OLIRef.*M;
        elseif M < 1 % Means OLI Ref is lower; so, multiply OLI Ref to scale it up
            OLI_ref_m = OLI_ref.*M;
            HyB_L8_m = HyBanded_OLIRef;
        end; clear M
        
     SSE_all_noscaling(ref) = sse(OLI_ref, HyBanded_OLIRef);
     SSE_all(ref) = sse(OLI_ref_m, HyB_L8_m); %clear OLI_ref_m HyB_L8_m
end; %clear ref OLI_ref HyBanded_OLIRef Hyper_WLenNM Hyper_Spectrum Hy_WL_Sand Hyper_Spectrum_all L8_allMean

%save('SBAF_L8_D1047_SCalSite2.mat'); % This file will have 3 variables in it