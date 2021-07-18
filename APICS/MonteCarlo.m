clear all
diff = randfixedsum(80,1,0.4525,-0.1,0.1);
load('L8andAPICS_Libya4')
load('Egyp1_OLI and Predicted.mat')
datenumber = date + datenum(num2str('2013-02-11'), 'yyyy-mm-dd');
%deciYear = decyear(datestr(datenumber));
deciYear = decimalYear;

%%
figure, plot(deciYear, diff*100, 'go',  'markers', 15, 'LineWidth', 2)
%RMSE = round(sqrt(mean((diff).^2)),5);
STD = round(std(diff),4);
tx = ['STD = ', num2str(100*STD), '%'];
text(2014, 13.5, tx, 'FontSize', 24,'Color','red')
ylim([-15 15])
%ylim([-15 15]), yticks([-5 -4 -3 -2 -1 0 1 2 3 4 5])
title('Pecentage Difference between Model and Satellite Observation')
ylabel('Percent Difference')
xlabel('Decimal Year')
grid minor; ax  = gca; ax.FontSize = 25; ax.GridColor = 'k';

%%
figure
for i = 1:1000
        for nlook = 1:10
            look_syn(nlook,i) = mean(randsample(diff, nlook, true));
        end
        plot(look_syn(:,i),'bo', 'MarkerSize', 10), hold on
end
%ylim([-0.14 0.14])

%%
Rand_Unc = std((look_syn),1,2);
mean_samples = mean((look_syn),2);
Rand_Unc_prc = round(100*Rand_Unc,2);
hold on
plot(std(look_syn,1,2), 'g--o', 'LineWidth', 3)
A(1:10) = Rand_Unc(5);
xlabel('Number of Observations')
ylabel('Differences and Uncertainty (Std. Dev.)')
title('Differences and Uncertainty')
tx = ['Green line tracks uncertainty which decreases with number of observations'];
text(1.5, 1.1*max(look_syn(1,:)),tx, 'FontSize', 24)

grid minor; ax  = gca; ax.FontSize = 25; ax.GridColor = 'k';
legend('Random Sample of Differences')

for i =1:9
    txt = ['Unc. = ', num2str(Rand_Unc_prc(i)),'%'];
    text(i+0.09, Rand_Unc(i)+0.003,txt,'Color','red','FontSize',16)
end

%% Percentage decrease in random error
Rand_Unc = std((look_syn),1,2);
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

%% How skewness affects the result
r = pearsrnd(2,1, -0.7, 3,10000,1);
figure, histogram(r)

%%
figure
for i = 1:10000%N(j)
        for nlook = 1:10
            look_syn2(nlook,i) = mean(randsample(r, nlook, true));
        end
        plot(look_syn2(:,i),'bo', 'MarkerSize', 10), hold on
end
%ylim([-0.14 0.14])
Rand_Unc = std((look_syn2),1,2);

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
