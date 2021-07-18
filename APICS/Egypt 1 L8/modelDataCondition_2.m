function [Data,DataTable,label,hyperionScaleFactor ] = modelDataCondition_2(satelliteNumber,...
    hyperionSpectrum, hyperionScaleFactor,hyperionWavelengths,sza,vza)
%UNTITLED2 Summary of this function goes here
% This function uses bander function. Define the satelliteNumber defined in
% the bander function.
% This function conditions the data and returns the rounded scale factor,
% BRDF coefficients.
% The function also requires the input of the gain factor and hyperion
% coincident pairs average spectrum.


% performing the banding to the hyperion data to the  satellite
[label,bandedReflectance] = bander(hyperionWavelengths,hyperionSpectrum,satelliteNumber);


% selecting the central wavelength of the satellite
[~,centralWavelength] = bander(hyperionWavelengths,hyperionWavelengths,...
    satelliteNumber);

% banding the scale factor to the satellites bands
[~,scaleFactor] = bander(hyperionWavelengths,hyperionScaleFactor, ...
    satelliteNumber);

% rounding the data to 5 decimal points 
scaleFactor = round(scaleFactor,6);
sqCentralWavelength = centralWavelength.^2;
cubeCentralWavelength = centralWavelength.^3;

% calculating the BRDF coefficient for each wavelength
BRDFsza = sza.a*exp(centralWavelength*sza.b) + sza.c*exp(centralWavelength*sza.d);
BRDFvzalin = vza.lin.a*sqCentralWavelength + vza.lin.b*centralWavelength + vza.lin.c;
BRDFvzaquad = vza.qua.a*sqCentralWavelength + vza.qua.b*centralWavelength + vza.qua.c;

Data(:,1) = BRDFsza;
Data(:,2) = BRDFvzalin;
Data(:,3) = BRDFvzaquad;
Data(:,5) = bandedReflectance;
Data(:,4) = scaleFactor;

% uncertainity table generation 
DataTable = table(BRDFsza',BRDFvzalin',BRDFvzaquad',bandedReflectance',scaleFactor');
DataTable.Properties.VariableNames = {'SZA','VZA_Lin','VZA_Qua','Band_Refl','Scale_Factor'};
end

