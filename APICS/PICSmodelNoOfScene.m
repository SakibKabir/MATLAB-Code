clear all;
load('Landsat8 Libya4NewAPICSFinal')

base='Z:\ImageDrive\OLI.TIRS\L8\P181\R040';
dates = dir(base);
dates([1 2])=[];

datenumber = date + datenum(num2str('2013-02-11'), 'yyyy-mm-dd');
deciYear = decyear(datestr(datenumber));

%%
L8_Red_band = Picsdata(:,9); % red band
Model_Red_band = Picsdata(:,21); 
date = Picsdata(:,16);
%%
figure, plot(date, L8_Red_band, 'ro', 'markers', 10)
hold on
plot(date, Model_Red_band, 'b*', 'markers', 10), ylim([0.32 0.36])
legend('Measured', 'Predicted')
title('Landsat 8 and Model Predicted TOA reflectance (Red Band)')
ylabel('TOA Reflectance')
xlabel('Day Since Launch')
grid minor; ax  = gca; ax.FontSize = 25; ax.GridColor = 'k';

%%
error = (L8_Red_band - Model_Red_band)./Model_Red_band;
figure, plot(date, error*100, 'go',  'markers', 15, 'LineWidth', 2)
ylim([-3 3])
title('Pecentage Difference L8 and Model Predicted TOA reflectance (Red Band)')
ylabel('Percentage Difference(%)')
xlabel('Day Since Launch')
grid minor; ax  = gca; ax.FontSize = 25; ax.GridColor = 'k';
RMSE = round(sqrt(mean((error).^2)),5);
STD = round(std(error),5);
tx = ['RMSE = ', num2str(100*RMSE), '%', ', STD = ', num2str(100*STD), '%'];
text(50, 1.8, tx, 'FontSize', 24)

%% Distribution of Error
figure,
histogram(error,30)
title('Error Distribution')
ylabel('Number of Observation')
xlabel('Error')
grid minor; ax  = gca; ax.FontSize = 25; ax.GridColor = 'k';
%%
meanError = mean(error);
absError = abs(error);
for nlook = 1:10
    for i = 1:30
        samples(nlook,i) = mean(randsample(error, 2, true));
    end
    sample_means = mean(samples); %clear samples
    sample_std = std(samples);
end

%%
figure,
for j=1:size(samples,2)
    plot(abs(samples(:,j)), '.', 'markers', 15)
    hold on
end
plot(abs(mean(samples,2)),'g--o')
%ylim([100*0.00292 100*0.00294]);
grid minor; ax  = gca; ax.FontSize = 35; ax.GridColor = 'k';
ylabel('Gain Estimate Error (%)')
xlabel('Number of Scene')

%% How RMSE change with number of data points
for i = 1:size(error,1)
   rmse(i,:) = rms(error(1:i,:));
   sd(i,:) = std(error(1:i,:));
end

