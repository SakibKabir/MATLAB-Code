clear all;
base = 'Z:\SpecialNeeds\Sakib\MATLAB\PICS Model Work';
base_dir = dir(base); base_dir([1,2]) = [];
file = fullfile(base_dir.folder, base_dir.name);
load(file) % loading the dataset
clearvars file base base_dir

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
grid minor; ax  = gca; ax.FontSize = 35; ax.GridColor = 'k';

%%
error = (L8_Red_band - Model_Red_band);
figure, plot(date, error*100, 'go', 'markers', 15)
ylim([-2 2])
title('Pecentage Difference L8 and Model Predicted TOA reflectance (Green Band)')
ylabel('Percentage Difference(%)')
xlabel('Day Since Launch')
grid minor; ax  = gca; ax.FontSize = 35; ax.GridColor = 'k';
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
        samples(i,nlook) = mean(randsample(error, nlook, true));
    end
    sample_means = mean(samples); %clear samples
    sample_std = std(samples);
end

%%
figure,
for j=1%:size(samples,1)
    plot(1:10, sample_std./abs(sample_means), '.', 'markers', 45)
    hold on
end

%ylim([100*0.00292 100*0.00294]);
grid minor; ax  = gca; ax.FontSize = 35; ax.GridColor = 'k';
ylabel('Gain Estimate Error (%)')
xlabel('Number of Scene')
