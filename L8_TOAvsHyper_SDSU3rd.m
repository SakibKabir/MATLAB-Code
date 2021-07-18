clear all
%%% Loading the Hyperspectral profiles and OLI Reflectance
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
        HyB_L8_m = HyB_L8.*M_vegi;
%         if M_vegi > 1 % Means OLI Ref is higer; so, multiply Hyper Ref to scale it up
%             OLI_ref_m = OLI_ref;
%             HyB_L8_m = HyB_L8.*M_vegi;
%         elseif M_vegi < 1 % Means OLI Ref is lower; so, multiply OLI Ref to scale it up
%             OLI_ref_m = OLI_ref.*M_vegi;
%             HyB_L8_m = HyB_L8;
%         end; clear HyB_L8
        % figure, plot(OLI_ref, '*r', 'MarkerSize', 18); hold on; plot(HyB_L8, 'pb', 'MarkerSize', 18)
        % figure, plot(OLI_ref_m, '*r', 'MarkerSize', 18); hold on; plot(HyB_L8_m, 'pb', 'MarkerSize', 18)
        if isnan(M_vegi)
            SSE(1, NoPro) = nan;
        else
            SSE(1, NoPro) = sse(OLI_ref, HyB_L8_m);
        end; clear M_vegi Hy_Spectrum_Vegi HyBanded_L8_Vegi OLI_ref_m HyB_L8_m
        
       if NoPro < 24
        %%% For Sand
        Hy_Spectrum_Sand = Hyper_Spectrum_all(size(Hyper_Spectrum, 1)+1:end, NoPro);
        %%% Banding Reference Sensor: OLI- Satnum = 24
        [Labels_L8_Sand, HyBanded_L8_sand] = bander(Hy_WL_Sand, Hy_Spectrum_Sand, 24); clear Hy_Spectrum_Sand Labels_L8_Sand
        HyB_L8s = HyBanded_L8_sand(2:5); clear HyBanded_L8_sand 
        M_sand = mean(OLI_ref./HyB_L8s);
        HyB_L8s_m = HyB_L8s.*M_sand;
%         if M_sand > 1
%             OLI_ref_m = OLI_ref;
%             HyB_L8_m = HyB_L8s.*M_sand;
%         elseif M_sand < 1
%             OLI_ref_m = OLI_ref.*M_sand;
%             HyB_L8_m = HyB_L8s;
%         end; clear M_sand HyB_L8s
%         
        % I Know I don't have any nan for sand
        SSE(2, NoPro) = sse(OLI_ref, HyB_L8s_m); 
       else 
       
       end 
       clear Labels_L8_Sand HyBanded_L8_sand OLI_ref_m HyB_L8_m
 end; %clear NoPro
 %%
SSE(SSE==0) = nan; 
%figure, histogram(SSE,10000); xlim([0 0.01])
%%
[r, c] = find(SSE < 0.01);
for sse = 1:size(r,1)
    r1 = r(sse); c1 = c(sse);
    s(sse,:) = SSE(r1, c1);
     if r1 == 1 % Vegi
            Hyper = Hyper_Spectrum_all(1:size(Hyper_Spectrum, 1), c1(1));
            Hyper_WL = Hyper_WLenNM; 
       elseif r1 == 2 % Sand
            Hyper = Hyper_Spectrum_all(size(Hyper_Spectrum, 1)+1:end, c1(1));
            Hyper_WL = Hy_WL_Sand;
     end %; clear r c 
%      Hyper_pro(sse,:) = {Hyper};
%      Hyper_WL_all(sse,:) = {Hyper_WL};
    % Banding 
    [labels_L8, HyBanded_L8] = bander(Hyper_WL, Hyper, 24); %clear labels_L8 
    HyBanded_OLIRef(sse,:) = HyBanded_L8(:, 2:5); %clear HyBanded_L8 

    %%% Banding Reference Sensor OLI- Satnum = 24
    [labels_L8, HyBanded_L8] = bander(Hyper_WL, Hyper, 24); clear labels_L8 
    HyBanded_OLIRef = HyBanded_L8(:, 2:5); clear HyBanded_L8
    %%% Banding Target Sensor Dove R(1047)- Satnum = 84
    [labels_1047, HyBanded_1047] = bander(Hyper_WL, Hyper, 84);  clear labels_1047 Hyper_WL Hyper 
    SBAF_L8_D1047_SCalSite(sse, :) = HyBanded_OLIRef./HyBanded_1047; clear HyBanded_1047
end
end
%%
figure
for sse = 1:size(s,1)
     plot(SBAF_L8_D1047_SCalSite(sse,:), '.','MarkerSize', 20)
     hold on
end

legend(['sse=', num2str(s(1))], ['sse=', num2str(s(2))])
figure, plot(s, SBAF_L8_D1047_SCalSite(:,4), '.')
%% Plotting SBAF vs SSE band by band
figure, 
for ii= 1:4
 figure(ii)
    plot(s*100, SBAF_L8_D1047_SCalSite(:,4), '.','MarkerSize', 30)
    grid on; grid minor; ax = gca; ax.FontSize = 30; ax.GridColor = 'k';
    ax.MinorGridColor = 'k'; % ax.GridAlpha = 0.8;
    title('SBAF vs Sum of Squared Error')
    xla
end


%%
% figure
% for i = 1:2%:size(Hyper_pro,1) 
%     H = cell2mat(Hyper_pro(i)); %cell2double(Hyper_pro()
%     WL = cell2mat(Hyper_WL_all(i));
%     plot(WL, H); xlim([0 1000])
%     hold on
% end
% m = mean(L8_allMean(1,:)./HyBanded_OLIRef(2,:));
% 
% figure, plot(m.*HyBanded_OLIRef(2,:), '.', 'MarkerSize', 20); ylim([0 1])
% hold on; plot(L8_allMean(1,:), 'r.', 'MarkerSize', 20);
% figure, plot(s, HyBanded_OLIRef(:,2), '.', 'MarkerSize', 20); ylim([0 1])
% grid on; grid minor; ax = gca; ax.FontSize = 30; ax.GridColor = 'k';
% ax.MinorGridColor = 'k'; % ax.GridAlpha = 0.8;