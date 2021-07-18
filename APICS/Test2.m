close all
clear all
load('L8andAPICS_Libya4')

datenumber = date + datenum(num2str('2013-02-11'), 'yyyy-mm-dd');
deciYear = decyear(datestr(datenumber));

%% Sample data Using Libya4 and Landsat 8 OLI
figure, plot(deciYear, L8_Red_band, 'ro', 'markers', 10)
hold on
plot(deciYear, Model_Red_band, 'b*', 'markers', 10), ylim([0.32 0.36])
legend('OLI measured', 'Model predicted')
title('Sample data Using Libya4 Model and Landsat 8 OLI (Red Band)')
ylabel('TOA Reflectance')
xlabel('Decimal Year')
grid minor; ax  = gca; ax.FontSize = 25; ax.GridColor = 'k';

%% Looking at the Percentage difference
diff = (Model_Red_band - L8_Red_band)./Model_Red_band;
figure, plot(deciYear, diff*100, 'go',  'markers', 15, 'LineWidth', 2)

RMSE = round(sqrt(mean((diff).^2)),5);
STD = round(std(diff),5);
tx = ['RMSE = ', num2str(100*RMSE), '%', ', STD = ', num2str(100*STD), '%'];
text(2013.3, 4.7, tx, 'FontSize', 24,'Color','red')
%ylim([-5 5]), yticks([-5 -4 -3 -2 -1 0 1 2 3 4 5])
title('Pecentage Difference between Model and L8 OLI (Red Band)')
ylabel('Percent Difference')
xlabel('Decimal Year')
grid minor; ax  = gca; ax.FontSize = 25; ax.GridColor = 'k';

%% Variability of the differences
%diff = (Model_Red_band - L8_Red_band)./Model_Red_band;
figure,
%subplot(1,2,1)
nbins = 20;
h = histogram(diff,nbins);
var_1sigma = round(std(diff), 4);
var_2sigma = round(2*std(diff),4);
var_3sigma = round(3*std(diff),4);
hold on
plot(mean(diff), 0 , 'r.', 'MarkerSize', 60)
e = errorbar(mean(diff), 0, var_1sigma, 'horizontal','LineWidth', 6);
%e.Marker = '*';
e.MarkerSize = 30;
e.Color = 'red';
e.CapSize = 15;
title('Distribution of Differences')
xlabel('Differences')
ylabel('Frequency')
grid minor; ax  = gca; ax.FontSize = 35; ax.GridColor = 'k';

tx = ['Standard Dev.: ', num2str(var_1sigma), ' (Red Bar)'];
text(-0.037, 24, tx, 'FontSize', 40, 'Color','red')
% tx = ['2SD: ', num2str(var_2sigma)];
% text(-0.035, 22.5, tx, 'FontSize', 20,'Color','blue')
% tx = ['3SD: ', num2str(var_3sigma)];
% text(-0.035, 21, tx, 'FontSize', 20,'Color','blue')

%% Looking at the Modeled estimates
% error = L8_Red_band-Model_Red_band; % Not using it anymore
% clear all
% load('L8andAPICS_Libya4')
diff = (Model_Red_band - L8_Red_band)./Model_Red_band;
for i = 1:10000
    for nlook = 1:10
        look(nlook,i) = mean(randsample(diff, nlook, true));
    end
    %plot(look(:,i),'bo', 'MarkerSize',10), hold on
end

%% Checking the effect of number of iteration
clear all
load('L8andAPICS_Libya4')
diff = (Model_Red_band - L8_Red_band)./Model_Red_band;
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

% var_1sigma =  round(std(std_Dev), 4);
% plot(mean(std_Dev), 'r.', 'MarkerSize', 15)
% e = errorbar(mean(std_Dev), 2*var_1sigma, 'vertical','LineWidth', 2);
% %e.Marker = '*';
% e.MarkerSize = 25;
% e.Color = 'red';
% e.CapSize = 15;

grid minor; ax  = gca; ax.FontSize = 25; ax.GridColor = 'k';
title('Standard Deviation vs Number of Iteration')
xlabel('Number of Iteration')
ylabel('Standard Deviation')

%% plotting histogram for 1 observation
j = 1;
figure, 
histogram(look(j,:))
title(['Distribution of Differences from 1 random sampling (10000 times)'])
xlabel('Estimated Differences')
ylabel('Frequency')
tx = ['Standard Dev.: ', num2str(round(std(look(j,:)),4)), ' (Red Bar)'];
text(-0.037, 2350, tx, 'FontSize', 40, 'Color','red')
hold on
plot(mean(look(j,:)), 0 , 'r.', 'MarkerSize', 60)
e = errorbar(mean(look(j,:)), 0, std(look(j,:)), 'horizontal','LineWidth', 6);
%e.Marker = '*';
e.MarkerSize = 30;
e.Color = 'red';
e.CapSize = 15;
grid minor; ax  = gca; ax.FontSize = 35; ax.GridColor = 'k';

%% plotting histogram for observation 2 to 5
figure,
for j =  2:5
    subplot(2,2,j-1)
    histogram(look(j,:))
    title([num2str(j),' Differences randomly selected and averaged (10000 times)'])
    xlabel('Estimated Differences')
    ylabel('Frequency')
    xlim([-0.04 0.04])
    grid minor; ax  = gca; ax.FontSize = 15; ax.GridColor = 'k';
    legend(['Stand. Dev.: ', num2str(round(100*std(look(j,:)),2)), '%'])
end

%% Plotting differences and Uncertainty
load('look')
figure,
for i = 1:10000
    plot(look(:,i),'bo', 'MarkerSize',10), hold on
end

%%
Rand_Unc = std((look),1,2);
mean_samples = mean((look),2);
Rand_Unc_prc = round(100*Rand_Unc,2);
hold on
plot(std(look,1,2), 'g--o', 'LineWidth', 3)
xlabel('Number of Observations')
ylabel('Differences and Uncertainty (Std. Dev.)')
title('Differences and Uncertainty')
tx = ['Green line tracks uncertainty which decreases with number of observations'];
text(1.5, 0.9*max(look(1,:)),tx, 'FontSize', 24)
grid minor; ax  = gca; ax.FontSize = 25; ax.GridColor = 'k';
legend('Random Sample of Differences')

for i =1:9
    txt = ['Unc. = ', num2str(Rand_Unc_prc(i)),'%'];
    text(i+0.09, Rand_Unc(i)+0.002,txt,'Color','red','FontSize',16)
end

%% Percentage decrease in random error
Rand_Unc = std((look),1,2);
SE_sx_hat = Rand_Unc;

%SE = std(diff);
for i = 1:size(Rand_Unc,1)
   % SE_all(i) = SE /sqrt(i);
    SE_sx_hat(i) = Rand_Unc(1)/sqrt(i);
end
%%
figure, plot(Rand_Unc,'b--o', 'LineWidth', 3), hold on
plot(SE_sx_hat,'g-', 'LineWidth', 3)
%ylim([0.03 0.06])
xlabel('Number of Observations')
ylabel('Standard Error and Uncertainty')
title('Uncertainty and Standard Error variation with Number of Observation')
grid minor; ax  = gca; ax.FontSize = 25; ax.GridColor = 'k';
legend('Uncertainty(from Monte Carlo)', 'Standard Errors')
