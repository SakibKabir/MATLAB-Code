
clear all
load('Egyp1_OLI and Predicted.mat')

L8_NIR_band = reflectance(:,5);
Model_NIR_band = predictedReflectance(:,5);

%% Sample data Using Libya4 and Landsat 8 OLI
figure, plot(decimalYear, L8_NIR_band, 'ro', 'markers', 10)
hold on
plot(decimalYear, Model_NIR_band, 'b*', 'markers', 10), ylim([0.52 0.62])
legend('OLI measured', 'Model predicted')
title('Sample data Using Libya4 Model and Landsat 8 OLI (Red Band)')
ylabel('TOA Reflectance')
xlabel('Decimal Year')
grid minor; ax  = gca; ax.FontSize = 25; ax.GridColor = 'k';

%% Looking at the Percentage difference
diff = (Model_NIR_band - L8_NIR_band)./Model_NIR_band;
figure, plot(decimalYear, diff*100, 'go',  'markers', 15, 'LineWidth', 2)
ylim([-10 10])
RMSE = round(sqrt(mean((diff).^2)),5);
STD = round(std(diff),5);
tx = ['RMSE = ', num2str(100*RMSE), '%', ', STD = ', num2str(100*STD), '%'];
text(2013.3, 4.7, tx, 'FontSize', 24,'Color','red')
%ylim([-5 5]), yticks([-5 -4 -3 -2 -1 0 1 2 3 4 5])
title('Pecentage Difference between Model and L8 OLI (Red Band)')
ylabel('Percent Difference')
xlabel('Decimal Year')
grid minor; ax  = gca; ax.FontSize = 25; ax.GridColor = 'k';

%% Checking the effect of number of iteration
L8_NIR_band = reflectance(:,5);
Model_NIR_band = predictedReflectance(:,5);
diff = (Model_NIR_band - L8_NIR_band)./Model_NIR_band;
iter = [10, 100, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000,...
        10000, 11000, 12000];

for a = 1:30
    for it = 1: length(iter)
        for i = 1:iter(it)
             look_num(i) = mean(randsample(diff, 2, true)); % 2 look
        end
        std_Dev(a, it) = std(look_num); clear look_num
    end
end

%% 
figure,
for a = 1:30
    plot(1:length(iter), std_Dev(a,:), 'o-', 'MarkerSize', 10), hold on
end

xticks(1: length(iter))
xticklabels({'10', '100', '1000', '2000', '3000', '4000', '5000', '6000',...
    '7000', '8000', '9000', '10000', '11000', '12000'})
var_1sigma =  round(std(std_Dev), 4);
plot(mean(std_Dev), 'r.', 'MarkerSize', 15)
e = errorbar(mean(std_Dev), var_1sigma, 'vertical','LineWidth', 2);
%e.Marker = '*';
e.MarkerSize = 25;
e.Color = 'red';
e.CapSize = 15;
grid minor; ax  = gca; ax.FontSize = 25; ax.GridColor = 'k';
title('Standard Deviation vs Number of Iteration')
xlabel('Number of Iteration')
ylabel('Standard Deviation')

%% Looking at the Modeled estimates
diff = (Model_NIR_band - L8_NIR_band)./Model_NIR_band;
figure
for i = 1:1000
    for nlook = 1:10
        look(nlook,i) = mean(randsample(diff, nlook, true));
    end
    plot(look(:,i),'bo', 'MarkerSize',10), hold on
end

%% Uncertainity Calculation
Rand_Unc = std((look),1,2);
mean_samples = mean((look),2);
Rand_Unc_prc = round(100*Rand_Unc,2);
Rand_Unc = std((look),1,2);

SE_sx_hat = Rand_Unc;
%SE = std(diff);
for i = 1:size(Rand_Unc,1)
   % SE_all(i) = SE /sqrt(i);
    SE_sx_hat(i) = Rand_Unc(1)/sqrt(i);
end
