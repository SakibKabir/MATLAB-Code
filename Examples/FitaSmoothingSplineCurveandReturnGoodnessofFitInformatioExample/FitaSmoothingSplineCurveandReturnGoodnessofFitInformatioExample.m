%% Fit a Smoothing Spline Curve and Return Goodness-of-Fit Information  
% Load some data and fit a smoothing spline curve through variables |month|
% and |pressure|, and return goodness of fit information and the output
% structure. Plot the fit and the residuals against the data.   

% Copyright 2015 The MathWorks, Inc.


%%  
load enso;
[curve, goodness, output] = fit(month,pressure,'smoothingspline');
plot(curve,month,pressure);
xlabel('Month');
ylabel('Pressure');  

%% 
% Plot the residuals against the x-data (|month|). 
plot( curve, month, pressure, 'residuals' )
xlabel( 'Month' )
ylabel( 'Residuals' )  

%% 
% Use the data in the |output| structure to plot the residuals against the
% y-data (|pressure|). 
plot( pressure, output.residuals, '.' )
xlabel( 'Pressure' )
ylabel( 'Residuals' )   

