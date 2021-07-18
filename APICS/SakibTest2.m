% A script to estimate reduction in error with multiple Landsat 8 looks
%  Make sure data file 'L8andAPICS_Libya4' is loaded
close all
clear all
load('L8andAPICS_Libya4')

%%
diff = L8_Red_band-Model_Red_band;
%figure
for j = 1:3%15%10000
%     for i = 1:30
%         look(1,i) = diff(round(117*rand+0.5));
%         look(2,i) = (diff(round(117*rand+0.5))+ diff(round(117*rand+0.5)))/2;
%         look(3,i) = (diff(round(117*rand+0.5))+ diff(round(117*rand+0.5))+diff(round(117*rand+0.5)))/3;
%         look(4,i) = (diff(round(117*rand+0.5))+ diff(round(117*rand+0.5))+diff(round(117*rand+0.5))+diff(round(117*rand+0.5)))/4;
%         look(5,i) = (diff(round(117*rand+0.5))+ diff(round(117*rand+0.5))+diff(round(117*rand+0.5))+diff(round(117*rand+0.5))+diff(round(117*rand+0.5)))/5;
%         look(6,i) = (diff(round(117*rand+0.5))+ diff(round(117*rand+0.5))+...
%         diff(round(117*rand+0.5))+diff(round(117*rand+0.5))+diff(round(117*rand+0.5))+diff(round(117*rand+0.5)))/6;
%         look(7,i) = (diff(round(117*rand+0.5))+ diff(round(117*rand+0.5))+diff(round(117*rand+0.5))+...
%         diff(round(117*rand+0.5))+diff(round(117*rand+0.5))+diff(round(117*rand+0.5))+diff(round(117*rand+0.5)))/7;
%         look(8,i) = (diff(round(117*rand+0.5))+ diff(round(117*rand+0.5))+diff(round(117*rand+0.5))+diff(round(117*rand+0.5))+...
%         diff(round(117*rand+0.5))+diff(round(117*rand+0.5))+diff(round(117*rand+0.5))+diff(round(117*rand+0.5)))/8;
%         %plot(abs(look(:,i)),'*'), hold on
%         %grid 
%         %pause
        
        for nlook = 1:5
            for i = 1:30
                look(nlook,i) = mean(randsample(diff, nlook, true));
            end
            
            %sample_means = mean(look); %clear samples
            %sample_std = std(look);
        end
        
    %end
    noLook{j,:} = look;
end

%%
figure,
for j=1:size(noLook,2)
    x = noLook{j,1};
    plot(abs(x(:,j)), '.', 'markers', 15)
    hold on
end
%plot(abs(mean(samples,2)),'g--o')
%ylim([100*0.00292 100*0.00294]);
grid minor; ax  = gca; ax.FontSize = 35; ax.GridColor = 'k';
ylabel('Gain Estimate Error (%)')
xlabel('Number of Scene')

%%
for j = 1:15%1000
    x = noLook{j,1};
    for ii = 1:nlook
        RMSE(j,ii) = sqrt(sum(abs(x(ii,:)))/30);
        SD(j,ii) = std(abs(x(ii,:))); 
    end
end

%% RMSE- which is the accuracy of the MODEL
figure
for j = 1:15%1000
    plot(RMSE(j,:), 'o'), hold on
end
plot(mean(RMSE), 'r--*')

%% SD- which is the accuracy of the MODEL
figure
for j = 1:15%1000
    plot(SD(j,:), 'o'), hold on
end
plot(mean(SD), 'r--*')

%%
%RMSE = round(sqrt(mean(look,2).^2),5);
plot(std(look,1,2), 'g--o')
%plot(abs(std(look,1,2)./mean(look,2)), 'r--o')
xlabel('Number of Observations')
ylabel('Absolute Difference and Uncertainty (Std. Dev.)')
title('Absolute Errors for 1-4 Observations and Uncertainty as a Function of Observation Number')
text(1.5,0.9*max(look(1,:)),'Green line tracks uncertainty which decreases with number of observations')
