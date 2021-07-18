%% To Compare with VIIRS 
clear all
load('c13_data.mat')

year_1 = c13_data(c13_data(:,11) <=365, :); nobs_y1 = size(year_1,1); % no of observation = 230
year_2 = c13_data(c13_data(:,11) <=730, :); nobs_y2 = size(year_2,1); % no of observation = 490
year_3 = c13_data(c13_data(:,11) <=1095, :);nobs_y3 = size(year_3,1); % no of observation = 756
year_4 = c13_data(c13_data(:,11) <=1460, :);nobs_y4 = size(year_4,1); % no of observation = 1016
year_5 = c13_data(c13_data(:,11) <=1825, :);nobs_y5 = size(year_5,1); % no of observation = 1254
year_6 = c13_data(c13_data(:,11) <=2190, :);nobs_y6 = size(year_6,1); % no of observation = 1494
year_7 = c13_data(c13_data(:,11) <=2555, :);nobs_y7 = size(year_7,1); % no of observation = 1730 % year 7 has 339 days

%%
[months_mean_1yr] = monthlyFunc(year_1);
[months_mean_2yr] = monthlyFunc(year_2);
[months_mean_3yr] = monthlyFunc(year_3);
[months_mean_4yr] = monthlyFunc(year_4);
[months_mean_5yr] = monthlyFunc(year_5);
[months_mean_6yr] = monthlyFunc(year_6);
[months_mean_7yr] = monthlyFunc(year_7);

%% Minimum Detectable Trend Calculation
[MDT_7, sigma_7, phi_7, acf_lags_bounds_7] = MDTCal(months_mean_7yr, 6.92);
[MDT_6, sigma_6, phi_6, acf_lags_bounds_6] = MDTCal(months_mean_6yr, 6);
[MDT_5, sigma5, phi_5, acf_lags_bounds_5] = MDTCal(months_mean_5yr, 5);
[MDT_4, sigma4, phi_4, acf_lags_bounds_4] = MDTCal(months_mean_4yr, 4);
[MDT_3, sigma3, phi_3, acf_lags_bounds_3] = MDTCal(months_mean_3yr, 3);
[MDT_2, sigma2, phi_2, acf_lags_bounds_2] = MDTCal(months_mean_2yr, 2);
[MDT_1, sigma1, phi_1, acf_lags_bounds_1] = MDTCal(months_mean_1yr, 1);

%% 
figure, plot(months_mean(:,2),'-o', 'MarkerSize', 10)

%%
acf_lags_bounds = acf_lags_bounds_7;
for b = 1:7
   acf = acf_lags_bounds{b}{1}; lags = acf_lags_bounds{b}{2};
   macf(:,b) = mean(acf,1);
   yt = months_mean_7yr(:,b);
   for i = 1:size(acf,1)-1
       yt_2 = acf(i+1).*yt(i);
       et(i,b) = yt(i) - yt_2;
   end
end

%%
figure, 
band_name ={'CA','Blue','Green','Red','NIR','SWIR1','SWIR2'};
for b = 1:7
    subplot(2,4,b),
    histogram(months_mean_7yr(:,b), 20)
    ylim([0 12])
    title(['ACF, ', band_name{b}, ' band']);
end

%% Looking at ACF
band_name ={'CA','Blue','Green','Red','NIR','SWIR1','SWIR2'};
figure,
acf_lags_bounds = acf_lags_bounds_7;
for b = 1:7
   acf = acf_lags_bounds{b}{1}; lags = acf_lags_bounds{b}{2};
   bounds = acf_lags_bounds{b}{3}; 
   subplot(2,4,b),
   stem(lags, acf, 'o', 'filled');ylim([-1 1])
   yline(bounds(1),'--r','LineWidth', 1.5); yline(bounds(2),'--r','LineWidth', 1.5);
   legend({'ACF of Cluster 13, L8','2-sigma confidence bounds'}, 'FontSize',12)
   title(['ACF, ', band_name{b}, ' band']);
   xlim([0 size(acf,1)+10])
   xlabel('Lags(months)'); ylabel('Sample Autocorrelation')
   hold on;  grid minor; ax  = gca; ax.FontSize = 12; ax.GridColor = 'k';
   ax.GridAlpha = 0.8; ax.MinorGridColor = 'k'; %ax.FontName = 'Times New Roman';
end
sgt = sgtitle('Sample Autocorrelation Function(ACF), near-daily observation(monthly averaged) for 1 year','Color','red');
sgt.FontSize = 20;


%%
%MDT_peryr_Raj = [1.83 1.66 1.34 1.70 2.00 2.60];
%%
band_name ={'Blue','Green','Red','NIR','SWIR1','SWIR2'};
X = categorical(band_name); X = reordercats(X, band_name);
figure,  plot(1:7, MDT_7(:,1:7), 'h--','MarkerSize', 10, 'LineWidth',1.5)
hold on, plot(1:7, MDT_6(:,1:7), '^--','MarkerSize', 10, 'LineWidth',1) 
hold on, plot(1:7, MDT_5(:,1:7), 's--','MarkerSize', 10, 'LineWidth',1.5) 
hold on, plot(1:7, MDT_4(:,1:7), 'd--','MarkerSize', 10, 'LineWidth',1)
hold on, plot(1:7, MDT_3(:,1:7), 'd--','MarkerSize', 10, 'LineWidth',1)
hold on, plot(1:7, MDT_2(:,1:7), 'p--','MarkerSize', 15, 'LineWidth',1)
hold on, plot(1:7, MDT_1(:,1:7), 's-','MarkerSize', 10, 'LineWidth',1.5)
yt = get(gca, 'XTick');

%text(yt, 1, labels, 'HorizontalAlignment','center', 'VerticalAlignment','bottom','FontSize',18,'Color','red')
hold on;  grid minor; ax  = gca; ax.FontSize = 24; ax.GridColor = 'k'; ax.GridAlpha = 0.8;
ax.MinorGridColor = 'k'; %ax.FontName = 'Times New Roman';
ylabel('Minimum Detectable Trend Per year (in %)')
title('Minimum Detectable Trend Comparison(L8 Cluster 13 and VIIRS Libya 4)')
legend({'VIIRS Libya 4 (9 Months) 1 obs. in 16 days', 'L8 C13 (9 Months) 1 obs. in 16 days',...
        'L8 C13 (9 Months) 1 obs. in 21 days',...
        'L8 C13 (9 Months) 1 obs. in 8 days','L8 C13 (9 Months) 1 obs. in 4 days',...
        'L8 C13 (9 Months) all obs.'}, 'Location','northwest' )

xticks([1 2 3 4 5 6])
xticklabels(band_name)

%% MDT change with sampling rate
MDT_srate = [MDT_1; MDT_2; MDT_3; MDT_4; MDT_5; MDT_6; MDT_7];

figure,
colors = {'b','c','g','r','m',[0.6667 0 0],'k' };
bands={'CA','Blue','Green','Red','NIR','SWIR1','SWIR2'};
mkr_mat=['o' '+' '*' '.' 'x' 's' 'd' '^' 'v' '>' '<' 'p' 'h' '<' 'p' 'h'];
xax = [1 2 3 4 5 6 6.92];
for b = 1:7
plot(xax, MDT_srate(:,b), 'marker', mkr_mat(b),'MarkerSize', 12, 'LineStyle',...
    '--','LineWidth',1,'color',colors{b}); hold on
end

hold on;  grid minor; ax  = gca; ax.FontSize = 24; ax.GridColor = 'k'; ax.GridAlpha = 0.8;
ax.MinorGridColor = 'k'; %ax.FontName = 'Times New Roman';
ylabel('Minimum Detectable Trend Per year (in %)')
xlabel('Length of the time series (years)'); 
title('L8 C13 Minimum Detectable Trend variation with Length of the Time Series')
legend(bands); %xlim([3 7])
yrs = {'1','2','3','4','5','6','6.92'};
xticks([1 2 3 4 5 6 6.92])
xticklabels(yrs)

%% Variability change with time
sigma = [sigma1; sigma2; sigma3; sigma4; sigma5; sigma_6; sigma_7];

figure,
colors = {'b','c','g','r','m',[0.6667 0 0],'k' };
bands={'CA','Blue','Green','Red','NIR','SWIR1','SWIR2'};
mkr_mat=['o' '+' '*' '.' 'x' 's' 'd' '^' 'v' '>' '<' 'p' 'h' '<' 'p' 'h'];
xax = [1 2 3 4 5 6 6.92];
for b = 1:7
plot(xax, sigma(:,b), 'marker', mkr_mat(b),'MarkerSize', 12, 'LineStyle',...
    '--','LineWidth',1,'color',colors{b}); hold on
end

hold on;  grid minor; ax  = gca; ax.FontSize = 24; ax.GridColor = 'k'; ax.GridAlpha = 0.8;
ax.MinorGridColor = 'k'; %ax.FontName = 'Times New Roman';
ylabel('Variability (in %)')
xlabel('Length of the time series (years)'); 
title('L8 C13 variability change with Length of the Time Series')
legend(bands); xlim([1 7]); ylim([0.2 2.2])
yrs = {'1','2','3','4','5','6','6.92'};
xticks([1 2 3 4 5 6 6.92])
xticklabels(yrs)

%% Autocorrelation change with time
phi = [phi_1; phi_2; phi_3; phi_4; phi_5; phi_6; phi_7];

figure,
colors = {'b','c','g','r','m',[0.6667 0 0],'k' };
bands={'CA','Blue','Green','Red','NIR','SWIR1','SWIR2'};
mkr_mat=['o' '+' '*' '.' 'x' 's' 'd' '^' 'v' '>' '<' 'p' 'h' '<' 'p' 'h'];
xax = [1 2 3 4 5 6 6.92];
for b = 1:7
plot(xax, phi(:,b), 'marker', mkr_mat(b),'MarkerSize', 12, 'LineStyle',...
    '--','LineWidth',1,'color',colors{b}); hold on
end

hold on;  grid minor; ax  = gca; ax.FontSize = 24; ax.GridColor = 'k'; ax.GridAlpha = 0.8;
ax.MinorGridColor = 'k'; %ax.FontName = 'Times New Roman';
ylabel('Autocorrelation')
xlabel('Length of the time series (years)'); 
title('L8 C13 Autocorrelation change with Length of the Time Series')
legend({'CA','Blue','Green','Red','NIR','SWIR1','SWIR2'}, 'Location','southeast' ); xlim([1 7]); ylim([-1 1])
yrs = {'1','2','3','4','5','6','6.92'};
xticks([1 2 3 4 5 6 6.92])
xticklabels(yrs)

%%
omega = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5] ;
for i = 1:size(omega,2)
for b = 1:7
phi_srt = sqrt((1+phi_7(:,b))/(1-phi_7(:,b)));
req_year(i,b) = (((2*sigma_7(:,b))/omega(i))*phi_srt).^(2/3);
end 
end

%%
colors = {'b','c','g','r','m',[0.6667 0 0],'k' };
bands={'CA','Blue','Green','Red','NIR','SWIR1','SWIR2'};
mkr_mat=['o' '+' '*' '.' 'x' 's' 'd' '^' 'v' '>' '<' 'p' 'h' '<' 'p' 'h'];

figure,
for b = 1:7
 plot(omega, req_year(:,b), 'marker', mkr_mat(b),'MarkerSize', 12, 'LineStyle',...
    '--','LineWidth',1,'color',colors{b}); hold on
end
hold on;  grid minor; ax  = gca; ax.FontSize = 24; ax.GridColor = 'k'; ax.GridAlpha = 0.8;
ax.MinorGridColor = 'k'; %ax.FontName = 'Times New Roman';
ylabel('Number of Years')
xlabel('Percetange of Trend per year'); 
title('Number of Years Required to Detect Certain Percent Trend')
legend({'CA','Blue','Green','Red','NIR','SWIR1','SWIR2'}, 'Location','northeast' ); 
xticks([0.1:0.1:1.5])
yticks([1:1:15])
%xlim([1 7]); ylim([-1 1])
% yrs = {'1','2','3','4','5','6','6.92'};
% xticks([1 2 3 4 5 6 6.92])
% xticklabels(yrs)

%% MDT change with record length
c13_data_temp_9months = c13_data(1:158,:); N = 9/12;
[months_mean_16_9m] = monthlyFunc(c13_data_temp_9months, 16);
[MDT_16_9m, sigma_16_9m, phi_16_9m] = MDTCal(months_mean_16_9m, N);

c13_data_temp_15months = c13_data(1:279,:); N = 15/12;
[months_mean_16_15m] = monthlyFunc(c13_data_temp_15months, 16);
[MDT_16_15m, sigma_16_15m, phi_16_15m] = MDTCal(months_mean_16_15m, N);

c13_data_temp_21months = c13_data(1:422,:); N = 21/12;
[months_mean_16_21m] = monthlyFunc(c13_data_temp_21months, 16);
[MDT_16_21m, sigma_16_21m, phi_16_21m] = MDTCal(months_mean_16_21m, N);

c13_data_temp_27months = c13_data(1:545,:); N = 27/12;
[months_mean_16_27m] = monthlyFunc(c13_data_temp_27months, 16);
[MDT_16_27m, sigma_16_27m, phi_16_27m] = MDTCal(months_mean_16_27m, N);

c13_data_temp_33months = c13_data(1:681,:); N = 27/12;
[months_mean_16_33m] = monthlyFunc(c13_data_temp_33months, 16);
[MDT_16_33m, sigma_16_33m, phi_16_33m] = MDTCal(months_mean_16_33m, N);

%%
band_name ={'Blue','Green','Red','NIR','SWIR1','SWIR2'};
X = categorical(band_name); X = reordercats(X, band_name);
figure, plot(1:6, MDT_peryr_Raj, 'o-', 'MarkerSize', 15, 'LineWidth',1.5)
hold on, plot(1:6, MDT_16_9m(:,2:7), 's--','MarkerSize', 15, 'LineWidth',1.5)
hold on, plot(1:6, MDT_16_15m(:,2:7), 'h--','MarkerSize', 15, 'LineWidth',1.5)
hold on, plot(1:6, MDT_16_21m(:,2:7), '^--','MarkerSize', 15, 'LineWidth',1.5) 
hold on, plot(1:6, MDT_16_27m(:,2:7), '^--','MarkerSize', 15, 'LineWidth',1.5)
hold on, plot(1:6, MDT_16_33m(:,2:7), 'p--','MarkerSize', 15, 'LineWidth',1.5) 

yt = get(gca, 'XTick');
%text(yt, 1, labels, 'HorizontalAlignment','center', 'VerticalAlignment','bottom','FontSize',18,'Color','red')
hold on;  grid minor; ax  = gca; ax.FontSize = 24; ax.GridColor = 'k'; ax.GridAlpha = 0.8;
ax.MinorGridColor = 'k'; %ax.FontName = 'Times New Roman';
ylabel('Minimum Detectable Trend Per year (in %)')
title('Minimum Detectable Trend Comparison(L8 Cluster 13 and Aqua MODIS Libya 4)')
legend({'VIIRS Libya 4 (9 Months) 1 obs. in 16 days', 'L8 C13 (9 Months) 1 obs. in 16 days',...
        'L8 C13 (9 Months) 1 obs. in 21 days',...
        'L8 C13 (9 Months) 1 obs. in 8 days','L8 C13 (9 Months) 1 obs. in 4 days',...
        'L8 C13 (9 Months) 1 obs. in 2 days'}, 'Location','northwest' )
   
xticklabels(band_name)

%% MDT change with time
MDT_all = [MDT_16_9m; MDT_16_15m; MDT_16_21m; MDT_16_27m; MDT_16_33m];
%%
figure, plot(MDT_all(:,1),'.-', 'MarkerSize', 15)

%%
clear all
load('c13_data.mat')

%%
i = 1;
c13_data_temp = c13_data;
for s_day = 1:30:max(c13_data_temp(:,11))
    months = c13_data_temp(c13_data_temp(:,11) >= s_day & c13_data_temp(:,11) <= s_day+29, :);
    if ~isempty(months)
      months_mean(i,:) = mean(months);
    else
      months_mean(i,:) = NaN(1, size(months,2));
    end
    i = i + 1;
end

[r,c]= find(isnan(months_mean));
months_mean(r,:) = [];

%% 1. Estimate sigma_n and phi from c13 daily observation
for b = 1:7
    sigma_n_monthly(:, b) = 100*(std(months_mean(:,b))/mean(months_mean(:,b)));
    [acf, lags, bounds] = autocorr(months_mean(:,b)); phi_monthly(:,b) = acf(2);
end

%% 2. Calculate time requirement in year (daily observation)
omega = [1 0.11 0.1 0.09 0.1 0.08 0.19]; % insert as percentage
N = max(c13_data_temp(:,11))/365;
for b=1:7
    %omega_daily(:,b) = (2*sigma_n_monthly(:,b)*(sqrt((1+phi_monthly(:,b))/(1-phi_monthly(:,b))))*(N^(-3/2)));
    req_year_monthly(:,b) = ((2*sigma_n_monthly(:,b)*(sqrt((1+phi_monthly(:,b))/(1-phi_monthly(:,b)))))/omega(b)).^(2/3);
end




