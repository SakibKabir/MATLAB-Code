% A script to estimate reduction in error with multiple Landsat 8 looks
%  Make sure data file 'L8andAPICS_Libya4' is loaded
close all
diff = L8_Red_band-Model_Red_band;
for i = 1:30
    look(1,i) = diff(round(117*rand+0.5));
    look(2,i) = (diff(round(117*rand+0.5))+ diff(round(117*rand+0.5)))/2;
    look(3,i) = (diff(round(117*rand+0.5))+ diff(round(117*rand+0.5))+diff(round(117*rand+0.5)))/3;
    look(4,i) = (diff(round(117*rand+0.5))+ diff(round(117*rand+0.5))+diff(round(117*rand+0.5))+diff(round(117*rand+0.5)))/4;
    plot(abs(look(:,i)),'*'), hold on
    grid 
    %pause
end
plot(std(look,1,2), 'g--o')
xlabel('Number of Observations')
ylabel('Absolute Difference and Uncertainty (Std. Dev.)')
title('Absolute Errors for 1-4 Observations and Uncertainty as a Function of Observation Number')
text(1.5,0.9*max(look(1,:)),'Green line tracks uncertainty which decreases with number of observations')
