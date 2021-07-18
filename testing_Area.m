%% Plotting OLI and OLI banded
figure,
m = mean(L8_allMean(6,:)./Hy_OLI_banded(6,:));
plot(L8_allMean(6,:), 'r*'); hold on; 
plot(Hy_OLI_banded(6,:), '.')

%%
clear all
[L8_allMean, L8_allSD, D_allMean, D_allSD] = OLI_DoveTOACCdata();

%% Checking the Sum of Squared Error
s = sort(SSE_all); 
ind_95 = round(0.95*length(s));
SNR_95prc = s(find(s) == ind_95);

ind_97 = round(0.97*length(s));
SNR_97prc = s(find(s) == ind_97);

figure,
h = histogram(SSE_all*100); % Histogram
N = max(h.Values); % Maximum bin count
SSE_all(SSE_all==0) = nan;
hold on
plot([SNR_95prc*100 SNR_95prc*100],[0 N/4],'g','LineWidth',3) % 95 percentage
plot([SNR_97prc*100 SNR_97prc*100],[0 N/4],'r','LineWidth',3) % 97 percentage
xlim([0 0.65])

xlabel('Sum of Squared Error in %')
ylabel('No. of Observation')
title('Sum of Squared Error Distribution');
hold on; grid on; grid minor; ax  = gca; ax.FontSize = 35; ax.GridColor = 'k';
%ax.FontName = 'Times New Roman'; %  ax.GridAlpha = 0.8; ax.MinorGridColor = 'k'; 
legend('Data', ['At 95% = ' num2str(round(SNR_95prc*100, 3)), ' percent'],...
       ['At 97% = ' num2str(round(SNR_97prc*100, 3)), ' percent'])

%% Checking the spatial uncertainty -- CV of Dove
UNC_NonU_b1 = D_allSD(:,1)./D_allMean(:,1);
UNC_NonU_b1_OLI = L8_allSD(:,1)./L8_allMean(:,1);
figure, histogram(UNC_NonU_b1_OLI)

UNC_NonU_b2 = D_allSD(:,2)./D_allMean(:,2);
figure, histogram(UNC_NonU_b2)

UNC_NonU_b3 = D_allSD(:,3)./D_allMean(:,3);
figure, histogram(UNC_NonU_b3)

UNC_NonU_b4 = D_allSD(:,4)./D_allMean(:,4);
figure, histogram(UNC_NonU_b4)