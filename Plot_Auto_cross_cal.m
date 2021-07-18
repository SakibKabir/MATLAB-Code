clear all
%% Plotting All bands 
band_name = {'BLUE','GREEN' ,'RED' ,'NIR' };
band_colors={'c','g','r','m'};
% Calling the function to extract data
[L8_allMean, L8_allSD, D_allMean, D_allSD] = OLI_DoveTOACCdata();
% load('SBAF_L8_D1047_SCalSite2.mat');
% D_allMean = D_allMean.*SBAF_L8_D1047_SCalSite;
%%
for band = 1:4
    figure(band), 
    errorbar(L8_allMean(:, band), D_allMean(:, band), D_allSD(:, band),'.','color', band_colors{band},'markers', 40)
    fitline = fit(L8_allMean(:, band), D_allMean(:, band), 'poly1');
    G_and_B = coeffvalues(fitline); Gain = round(G_and_B(1,1),4); Bias = round(G_and_B(1,2),4);
    hold on; plot([0:0.0001:0.7], Gain*[0:0.0001:0.7]+Bias, 'b','LineWidth', 1); xlim([0 0.7]); ylim([0 0.7])
    Gain_all(band) = Gain; Bias_all(band) = Bias; ci = confint(fitline); ci = round(ci, 4);
    Gain_LB = ci(1,1); Bias_LB = ci(1,2); Gain_UB = ci(2,1); Bias_UB = ci(2,2);
    
    mdl = fitlm(L8_allMean(:,band), D_allMean(:,band));
    R_squared = round(mdl.Rsquared.Ordinary, 4); RMSEprc = round(mdl.RMSE*100,4);  
    NoofObs = mdl.NumObservations ;
    hold on; plot([0 0.7], [0 0.7], '--k', 'LineWidth', 0.5)
    
    equation = [num2str(Gain) '*x + ' num2str(Bias)]; 
    tx= strcat('y =  ', '  ', equation, ' (Fitted Blue Line)');
    text(0.01, 0.68, tx, 'FontSize', 24);
    tx = ['No. of Observation = ', num2str(NoofObs), '; R^{2} = ', num2str(R_squared),...
         ', RMSE = ', num2str(RMSEprc), '%'];
    text(0.01, 0.65, tx, 'FontSize', 24)
    
    tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(Gain),')',...
         ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
    text(0.01, 0.62, tx, 'FontSize', 24);
     
    tx = strcat('Bias', '(', num2str(Bias), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
    text(0.01   , 0.59, tx, 'FontSize', 24); 
    
%     tx = strcat('Blue line: Lower bound;', ' Green line: Upper bound');
%     text(0.02, 0.53, tx, 'FontSize', 18);

            
%     hold on;
%     %figure,
%     plot([0, 3], Gain_LB*[0, 3]+Bias_LB, '-b')
    
%     hold on;
%     plot([min(L8_allMean(:,band)), max(L8_allMean(:,band))], Gain_UB*[min(L8_allMean(:,band)), max(L8_allMean(:,band))]+Bias_UB,'-g')
%      
    %         plot(L8_band2, D_band1,'.',...
%         [min(L8_band2),max(L8_band2)],ci(1,1)*[min(L8_band2),max(L8_band2)]+ci(1,2),'-b',...
%         [min(L8_band2),max(L8_band2)],ci(2,1)*[min(L8_band2),max(L8_band2)]+ci(2,2),'-g')

%     plot(L8_allMean(:,band), D_allMean(:,band),'.',...
%         [min(L8_allMean(:,band)), max(L8_allMean(:,band))], ci(1,1)*[min(L8_allMean(:,band)), max(L8_allMean(:,band))]+ci(1,2),...
%         [min(L8_allMean(:,band)), max(L8_allMean(:,band))], ci(2,1)*[min(L8_allMean(:,band)), max(L8_allMean(:,band))]+ci(2,2))
    
    title(strcat('Mean TOA Reflectance Comparison', {', '} , band_name{band},' Band'));
    xlabel('Mean TOA Reflectance of L8 ROI')
    ylabel('Mean TOA Reflectance of D1047 ROI')
    hold on; grid on; grid minor; ax  = gca; ax.FontSize = 35; ax.GridColor = 'k';
    %ax.FontName = 'Times New Roman'; %  ax.GridAlpha = 0.8; ax.MinorGridColor = 'k'; 
end; clear band Gain Bias G_and_B fitline tx tx


