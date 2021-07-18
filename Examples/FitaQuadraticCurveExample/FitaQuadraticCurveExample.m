%% Fit a Quadratic Curve  
% Load some data, fit a quadratic curve to variables |cdate| and |pop|,
% and plot the fit and data.   

% Copyright 2015 The MathWorks, Inc.


%%  
load census;
f=fit(cdate,pop,'poly2')
plot(f,cdate,pop) 

%%
% For a list of library model names, see |fitType|.   

