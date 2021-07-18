clear all
[Hyper_WLenNM, Hyper_Spectrum, Hy_WL_Sand,Hyper_Spectrum_all] = HyperData();
[L8_allMean, D_allMean, D_allSD] = OLI_DoveTOACCdata();

%% L8_allMean = OLI_DoveTOACCdata();
for ref = 1 %: size(L8_allMean) % all the reflectances inside the ROI
    OLI_ref = L8_allMean(ref,:);
for NoPro = 9:size(Hyper_Spectrum_all, 2)
        %% For Vegitation
        Hy_Spectrum_Vegi = Hyper_Spectrum_all(1:size(Hyper_Spectrum,1), NoPro);
        %Hy_WL_Vegi = Hyper_WLenNM;
        %%% Banding Reference Sensor: OLI- Satnum = 24
        [Labels_L8_Vegi, HyBanded_L8_Vegi] = bander(Hyper_WLenNM, Hy_Spectrum_Vegi, 24);
        HyB_L8 = HyBanded_L8_Vegi(2:5); %clear Labels_L8_Vegi HyBanded_L8_Vegi
        
        figure, plot(OLI_ref, '*r', 'MarkerSize', 18); hold on; plot(HyB_L8,'pb', 'MarkerSize', 18)
        S_factor = mean(OLI_ref)/mean(HyB_L8);
       % figure, plot(OLI_ref, '*r', 'MarkerSize', 18); hold on; plot(HyB_L8.*S_factor,'pb', 'MarkerSize', 18)
        
        %HyB_L8 = HyB_L8.*S_factor;
        
        Diff = OLI_ref - HyB_L8;
        S_OLIRef = [OLI_ref(2)-OLI_ref(1), OLI_ref(3)-OLI_ref(2), OLI_ref(4)-OLI_ref(3)];
        S_HyB_L8 = [HyB_L8(2)-HyB_L8(1), HyB_L8(3)-HyB_L8(2), HyB_L8(4)-HyB_L8(3)];

        if (Diff(1)>0 &&  Diff(2)>0 && Diff(3)>0 && Diff(4)>0) || (Diff(1)<=0 &&  Diff(2)<=0 && Diff(3)<=0 && Diff(4)<= 0)
            if (S_OLIRef(1) >0 && S_HyB_L8(1)> 0) && (S_OLIRef(2) >0 && S_HyB_L8(2)> 0)&& (S_OLIRef(3) >0 && S_HyB_L8(3)> 0)
                S_factor = mean(OLI_ref)/mean(HyB_L8);
                
                Slope_Diff(1, NoPro) = sse(S_OLIRef, S_HyB_L8);
                % Slope_ratio(1, NoPro) = S_OLIRef./S_HyB_L8;
            elseif (S_OLIRef(1) <0 && S_HyB_L8(1)< 0) && (S_OLIRef(2) <0 && S_HyB_L8(2)< 0)&& (S_OLIRef(3) <0 && S_HyB_L8(3)< 0)
                Slope_Diff(1, NoPro) = sse(S_HyB_L8, S_OLIRef);
                Slope_ratio(1, NoPro) = S_HyB_L8/S_OLIRef;
            else
                Slope_Diff(1, NoPro) = nan;
                Slope_ratio(1, NoPro) = nan;
            end
            
        else
                Slope_Diff(1, NoPro) = nan;
                Slope_ratio(1, NoPro) = nan;
        end 
        clear Hy_Spectrum_Vegi Labels_L8_Vegi HyBanded_L8_Vegi 
        
        %% For Sand
        Hy_Spectrum_Sand = Hyper_Spectrum_all(size(Hyper_Spectrum,1)+1:end, NoPro);
        %%% Banding Reference Sensor: OLI- Satnum = 24
        [Labels_L8_Sand, HyBanded_L8_sand] = bander(Hy_WL_Sand, Hy_Spectrum_Sand, 24);
        HyB_L8S = HyBanded_L8_sand(2:5); 
        
        S_factors = mean(OLI_ref)/mean(HyB_L8S);
        HyB_L8S = HyB_L8S.*S_factors;
        % figure, plot(OLI_ref, '*r', 'MarkerSize', 18); hold on; plot(HyB_L8S,'pb', 'MarkerSize', 18)
        
        Diff = OLI_ref - HyB_L8S;
        S_OLIRef = [OLI_ref(2)-OLI_ref(1), OLI_ref(3)-OLI_ref(2), OLI_ref(4)-OLI_ref(3)];
        S_HyB_L8 = [HyB_L8S(2)-HyB_L8S(1), HyB_L8S(3)-HyB_L8S(2), HyB_L8S(4)-HyB_L8S(3)];
        
        if (Diff(1)>0 &&  Diff(2)>0 && Diff(3)>0 && Diff(4)>0) || (Diff(1)<=0 &&  Diff(2)<=0 && Diff(3)<=0 && Diff(4)<= 0)
            if (S_OLIRef(1) >0 && S_HyB_L8(1)> 0) && (S_OLIRef(2) >0 && S_HyB_L8(2)> 0)&& (S_OLIRef(3) >0 && S_HyB_L8(3)> 0)
                Slope_Diff(2, NoPro) = sse(S_OLIRef - S_HyB_L8);
                Slope_ratio(2, NoPro) = S_OLIRef/S_HyB_L8;
            elseif (S_OLIRef(1) <0 && S_HyB_L8(1)< 0) && (S_OLIRef(2) <0 && S_HyB_L8(2)< 0)&& (S_OLIRef(3) <0 && S_HyB_L8(3)< 0)
                Slope_Diff(2, NoPro) = sse(S_HyB_L8 - S_OLIRef);
                Slope_ratio(2, NoPro) = S_HyB_L8/S_OLIRef;
            else
                Slope_Diff(2, NoPro) = nan;
                Slope_ratio(2, NoPro) = nan;
            end
            
        else
                Slope_Diff(2, NoPro) = nan;
                Slope_ratio(2, NoPro) = nan;
        end    
        clear Hy_Spectrum_Sand Labels_L8_Sand HyBanded_L8_sand 
end; clear NoPro       
    
%%
Slope_ratio2 = abs(1-Slope_ratio);
[r1,c1] = find(Slope_ratio2 == min(min(Slope_ratio2)));
[r,c] = find(Slope_Diff == min(min(Slope_Diff)));

%%
    if isempty(r1)
        SBAF_L8_D1047_SCalSite(ref, :) = nan(1,4);
    else    
        if r1(1) == 1 % Vegi
            Hyper = Hyper_Spectrum_all(1:size(Hyper_Spectrum, 1), c1(1));
            Hyper_WL = Hyper_WLenNM;

        elseif r1(1) == 2 % Sand
            Hyper = Hyper_Spectrum_all(size(Hyper_Spectrum, 1)+1:end, c1(1));
            Hyper_WL = Hy_WL_Sand;
        end

        %%% Banding Reference Sensor OLI- Satnum = 24
        [labels_L8, HyBanded_L8] = bander(Hyper_WL, Hyper, 24);
        %%% Banding Target Sensor Dove R(1047)- Satnum = 84
        [labels_1047, HyBanded_1047] = bander(Hyper_WL, Hyper, 84);
        SBAF_L8_D1047_SCalSite(ref, :) = HyBanded_L8(:, 2:5)./HyBanded_1047;

    end 
end   
    
%%
figure, plot(OLI_ref, '*' )
hold on; plot(HyB_L8S, 'p')
% A = [1 2 1 2];
% B = [2.1 3.1 2.1 1]; 
% x1 = A - B; 
% figure, plot(B, '*'); xlim([0 5]); ylim([0 5])
% hold on; plot(A, 'o'); xlim([0 5]); ylim([0 5])
% [p_slope1, s1] = polyfit([1:4], B, 3);
% f1 = polyval(p_slope1, [1:4]);
% hold on;
% plot([1:4], f1)
% 
% [p_slope2, s2] = polyfit([1:4], A, 3);
% f2 = polyval(p_slope2, [1:4]);
% hold on;
% plot([1:4], f2)
% 
% 
% p_slope1(1,1:3) - p_slope2(1,1:3)
% x = fit([1:4]',B','poly1');