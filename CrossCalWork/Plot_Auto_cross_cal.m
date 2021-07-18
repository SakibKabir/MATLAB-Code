%% Plotting Automatic ROI 
band_colors={'c','g','r','m'};
band_name = {'BLUE', 'GREEN', 'RED', 'NIR'};
band = 2;

figure, plot(L8_meanTOAref_ROI_5B, S2A_meanTOAref_ROI_5B, '.', 'color', band_colors{band},'markers', 30)
%xlim([0.2 0.4])
%ylim([0.2 0.4])
% hold on
% plot([0 0.4], [0 0.4], 'k')
  
% Radiance
title(strcat('Mean TOA Reflectance Comparison', {', '} , band_name{band},' Band'));
xlabel('Mean TOA Reflectance of L8 (W/sr/m^2/{\mum})')
ylabel('Mean TOA Reflectance of D1047 (W/sr/m^2/{\mum})')

hold on
grid on
grid minor
ax  = gca;
ax.FontSize = 25;
ax.GridColor = 'k';
%   ax.GridAlpha = 0.8;
ax.MinorGridColor = 'k';
%ax.FontName = 'Times New Roman';

%%
p_slope = polyfit(L8_meanTOAref_ROI_5G, Dove_meanTOAref_ROI_5G, 1);
p_slope = round(p_slope, 5);
fitline_ind = fit(L8_meanTOAref_ROI_5G(:), Dove_meanTOAref_ROI_5G(:), 'poly1');
idx_d = find(~isnan(Dove_meanTOAref_ROI_5G(:)));
Dove_green_temp = Dove_meanTOAref_ROI_5G(:);
Dove_green = Dove_green_temp(idx_d);

%idx = find(~isnan(L8_meanTOAref_ROI_5G(:)));
OLI_green_temp = L8_meanTOAref_ROI_5G(:);
OLI_green = OLI_green_temp(idx_d);

fitline_ind = fit(OLI_green, Dove_green, 'poly1');

