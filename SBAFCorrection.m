clc 
clear 
close all
%%
% Libya 4 data Selection only
selcol = [6,7,8,19,21];

%{
%}
% performing the banding of Hyperion over Terra Modis Bands
% Terra Modis contains 7 Bands- CA, Blue, Green, Red NIR, SWIR1, SWIR2
modisBands = 7;

% Loading the coincident pairs data
load('C:\Users\bipin.raut\Documents\MATLAB\Modis_Hyperion_Libya_4_coincident.mat')

% loading the hyperion Data 
load('Z:\SpecialNeeds\BIPIN RAUT\ReflectanceData\Libya 4\Markot Libya ROI\hyperion_over_Libya_4.mat')
%load('C:\Users\bipin.raut\Documents\MATLAB\mat files\hyperion_over_Libya_4.mat')

% loading the modis Reflectance for the 
[modisReflectance,solarAngle,sensorAngle] = modis_txt('Libya4_ONLY_MODIS.dat');

%Drift correction of each band of Hyperion Using Xin Jing Data
driftPerYear = xlsread('C:\Users\bipin.raut\Documents\MATLAB\excel files\driftPerYear');
% converting into data set  for the ease of access
driftPerYear = dataset({driftPerYear, 'Wavelength','driftPerYear','slope','bias'});

%using the ROI provided by the Markot
%load('Z:\SpecialNeeds\BIPIN RAUT\ReflectanceData\Libya 4\Markot Libya ROI\hyperion_over_Libya_4.mat')

%###############################
%summaryAcquisitionDate(3,:) = [];
summaryAcquisitionDate = summaryAcquisitionDate(selcol,:);

% Extracting the solar zenith and solar azimuth  angles for all the co-incident pairs of
% Hyperion
hyperionSolarZenithCoincident = sunZenith(summaryAcquisitionDate.Hyperion_Row);
hyperionSolarAzimuthCoincident = sunAzimuth(summaryAcquisitionDate.Hyperion_Row);
hyperionSensorZenithCoincident = sensorLookAngle(summaryAcquisitionDate.Hyperion_Row);

% separating the wavelengths from the meanReflectance variable
% first column consist of the wavelength
hyperionWavelengths = meanReflectance(1,:);

% seperating the reflectance of hyperion from variable 'meanReflectance'
hyperionReflectance = meanReflectance(2:end,:);

% selecting the Hyperion Reflectance of the co-incident pair
hyperionReflectanceCoincident = hyperionReflectance(summaryAcquisitionDate.Hyperion_Row,:);


% plotting the hyperion selected acquistions
labels = num2str((1:length(summaryAcquisitionDate.Hyperion_Row))');
figure 
for sel = 1:length(summaryAcquisitionDate.Hyperion_Row)
    plot(hyperionWavelengths,hyperionReflectanceCoincident(sel,:));
    title('Hyperion Spectrum of Coincident pairs');
    xlabel('Wavelengths');
    ylabel('Reflectance');
    text(1600,hyperionReflectanceCoincident(sel,118),labels(sel));
    hold on
    grid on
    ax = gca;
    ax.FontSize = 20;
end

averageHyperionReflectanceCoincident = mean(hyperionReflectanceCoincident);
stdHyperionReflectanceCoincident = (std(hyperionReflectanceCoincident)).*100;

% Initalizing the BandedHyperionReflectance
% Initalization uses the count of acquisition which is row length of
% "summaryAcquisitionDate and "modisBands" for column
bandedHyperionReflectance = zeros(length(summaryAcquisitionDate.Hyperion_Row)...
    ,modisBands);

% correcting for the drift in the hyperion
hyperionCoincidentPairsDsl = summaryAcquisitionDate.HyperionDSL;

[coincidentPairNumber,hyperionBands] = size(hyperionReflectanceCoincident);

% repeating the matrix to the length of the reflectacnce to make
% calculation straight forward
hyperionConcidentPairsDslRepeat = repmat(hyperionCoincidentPairsDsl,1,...
    hyperionBands);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% applying the drift correction formula on th reflectance

% FORMULA:
%           rho_drift_corrected = rho - ((drift%peryear)*DSL/(365*100))

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% hyperionReflectanceCoincidentDriftCorrected = hyperionReflectanceCoincident - ...
%     (repmat(transpose(driftPerYear.driftPerYear),coincidentPairNumber,1) ...
%     .*hyperionConcidentPairsDslRepeat)/(365*100);
% 
% %plotting the drift corrected and not corrected hyperion reflectance in the same plot
% color = {'r','g','b','m','c'};
% figure 
% for sel = 1:length(summaryAcquisitionDate.Hyperion_Row)
% %     nonCorrected = plot(hyperionWavelengths,hyperionReflectanceCoincident(sel,:));
% %     nonCorrected.Color = color{sel};
% %     nonCorrected.Marker = 'o';
%     title('Drift Corrected Reflectance ');
%     xlabel('Wavelengths');
%     ylabel('Reflectance');
%     text(1600,hyperionReflectanceCoincident(sel,118),labels(sel));
%     hold on
%     corrected = plot(hyperionWavelengths,hyperionReflectanceCoincidentDriftCorrected(sel,:));
%     %corrected.Color = color{sel};
% %     corrected.Marker = '*';
% end
% 
% %{
% Applying the calibration Correction Coefficients Validation 
% Author : Xin Jing, Larry Leigh and Dennis Helder
% Paper : Lifetime Absolute Calibration of the EO-1 Hyperion Sensor and its
%         validation
% 
% The correction coefficeints : gain and bias
% 
% FORMULA :
%                 rho = rho_drift_corrected *slope + bias
% 
% %}
% hyperionReflectanceCoincidentDriftAndGainCorrected = ...
%     repmat(transpose(driftPerYear.slope),coincidentPairNumber,1)...
%     .*hyperionReflectanceCoincidentDriftCorrected ...
%     + repmat(transpose(driftPerYear.bias),coincidentPairNumber,1);
% 
% % plotting non corrected drift corrected gain corrected and non corrected
% % data in the same plot 
% color = {'r','g','b','m','c'};
% figure 
% for sel = 1:length(summaryAcquisitionDate.Hyperion_Row)
% %     nonCorrected = plot(hyperionWavelengths,hyperionReflectanceCoincident(sel,:));
% %     nonCorrected.Color = color{sel};
% %     nonCorrected.Marker = 'o';
% %     nonCorrected.LineStyle = '--';
%     title('Gain and Bias correction after Drift correction ');
%     xlabel('Wavelengths');
%     ylabel('Reflectance');
%     text(1600,hyperionReflectanceCoincident(sel,118),labels(sel));
%     hold on
% %     driftCorrected = plot(hyperionWavelengths,hyperionReflectanceCoincidentDriftCorrected(sel,:));
% %     driftCorrected.Color = color{sel};
% %     driftCorrected.Marker = '*';
% %     driftCorrected.LineStyle = '-.';
%     hold on
%     gainCorrected = plot(hyperionWavelengths,hyperionReflectanceCoincidentDriftAndGainCorrected(sel,:));
%    % gainCorrected.Color = color{sel};
%     %gainCorrected.Marker = '^';
%     %gainCorrected.LineStyle = ':';
% end
% 


% using the bander function to calculate the reflectance of the Hyperion
% using TerraModis RSR
for sel = 1:length(summaryAcquisitionDate.Hyperion_Row)
    [labels,Banded] = bander(hyperionWavelengths,hyperionReflectanceCoincident(sel,:),26);
    bandedHyperionReflectance(sel,:) = Banded;
end



% clearing the memory by removing the unused variables
%clear Banded hyperionReflectance hyperionReflectanceCoincident hyperionWavelengths ...
    %sunAzimuth sunElevation sunZenith ;



% Extracting the solar zenith and solar azimuth  angles for all the co-incident pairs of
% Hyperion
modisSolarCoincident = solarAngle(summaryAcquisitionDate.Modis_Row,:);
modisSensorCoincident = sensorAngle(summaryAcquisitionDate.Modis_Row,:);

% selecting the modis reflectance of the co-incident pair
modisReflectanceCoincident = modisReflectance(summaryAcquisitionDate.Modis_Row,:);

% formatting the Data to modis for comparision
modisReflectanceCoincidentArranged = [modisReflectanceCoincident(:,3),...
    modisReflectanceCoincident(:,4), modisReflectanceCoincident(:,1),...
    modisReflectanceCoincident(:,2),modisReflectanceCoincident(:,6),...
    modisReflectanceCoincident(:,7)];

% calculating the correction factor. The correction factor is given by the
% division of the modis divided by the hyperion reflectance
scaleFactor = modisReflectanceCoincidentArranged ./ bandedHyperionReflectance(:,2:end);
scaleFactorMean = mean(scaleFactor);
%scaleFactorMean = [0.9794,1.0092,1.0163,0.9844,0.9923,0.9681];

scaleFactorStd = std(scaleFactor);

%center wavelength of the Terra Modis
[~,modisCenterWavelength] = bander(hyperionWavelengths,hyperionWavelengths,26);
modisCenterWavelength = modisCenterWavelength(2:end);
% center wavelength of the Costal Aerosol
CAwavelength = 441.5;

% extrapolating the scaling factor to the costal band center wavelength
CAMeanScalingFactor = interp1(modisCenterWavelength,scaleFactorMean...
    ,CAwavelength,'Linear','extrap');

% extrapolation of the standard deviation 
CAStdScalingFactor = interp1(modisCenterWavelength,scaleFactorStd...
    ,CAwavelength,'Linear','extrap');

% adding the extrapolated data to the mean scale factor
ScaleFactorMeanWithoutCA = scaleFactorMean;
scaleFactorMean = [CAMeanScalingFactor,scaleFactorMean];

%adding the extrapolates standard deviation to the scaleFactoestd
scaleFactorStd = [CAStdScalingFactor,scaleFactorStd];

%adding not available data to the first coulmn of the scaleFactor
CAEmptyColumn = nan(sel,1);
scaleFactor = [CAEmptyColumn,scaleFactor];

% selecting coincident days 
coincidentDates = cellstr(datestr(summaryAcquisitionDate.Modis_Acquisition_Libya_1,'yyyy-mm-dd'));

% addining all the information in one table
NormalizationInfo = [scaleFactor;scaleFactorMean;scaleFactorStd];

% creating a data frame for the scale fator and its mean and standard
% deviation
scaleFactorDataSet = dataset({NormalizationInfo 'CA','Blue', 'Green','RED', 'NIR',...
    'SWIR1','SWIR2'});
scaleFactorDataSet.Properties.ObsNames = [coincidentDates;{'Mean';'STD'}];

% angle information for each co-incident date 
coincidentPairsAngleSummary = [hyperionSensorZenithCoincident,...
    hyperionSolarZenithCoincident, hyperionSolarAzimuthCoincident ...
    modisSensorCoincident(:,1), modisSensorCoincident(:,2),...
    modisSolarCoincident(:,1),modisSolarCoincident(:,2)];

% angle mean information  for each co-incident date.
meanAngle = [mean(hyperionSensorZenithCoincident),mean(hyperionSolarZenithCoincident),...
    mean(hyperionSolarAzimuthCoincident),mean(modisSensorCoincident(:,1)),...
    mean(modisSensorCoincident(:,2)),mean(modisSolarCoincident(:,1)),mean(modisSolarCoincident(:,2))];

% concatenating the mean information of the angles
coincidentPairsAngleSummary = [coincidentPairsAngleSummary;meanAngle];

% Finally interpolating the scale Factor to all the 



% creating the a dataFrame for the angle summary for the co-incident pairs 
coincidentPairsAngleSummary = dataset({coincidentPairsAngleSummary ...
    'Hyperion_Sensor' ,'Hyperion_Solar_Zenith', 'Hyperion_Solar_Azimuth', ...
    'Modis_Sensor_Zenith', 'Modis_Sensor_Azimuth', 'Modis_Solar_zenith',...
    'Modis_Solar_Azimuth'});
coincidentPairsAngleSummary.Properties.ObsNames = [coincidentDates;{'MEAN'}];

% interpolating the scaling Factor 
gainFactor =  interp1(modisCenterWavelength,ScaleFactorMeanWithoutCA, ...
    hyperionWavelengths,'linear','extrap');


%{
    Interpolating between the bands with different techniques to get gains
   %}
    
 % selecting the hyperion Bands to assign the gains
CARegion = find(hyperionWavelengths > 400 & hyperionWavelengths <= 500);
blueGreenRegion = find(hyperionWavelengths >  500 & hyperionWavelengths <= 590);
RedRegion = find(hyperionWavelengths > 590 & hyperionWavelengths <= 835);
nirRegion = find(hyperionWavelengths > 835 & hyperionWavelengths <= 960);
swir1Region = find(hyperionWavelengths >  960 & hyperionWavelengths <=1850);
swir2Region = find(hyperionWavelengths > 1850 & hyperionWavelengths <= 2395);

% BandGains for each variables
% initializing the varibale Bandgains with zeros
bandGain = zeros(1,hyperionBands);
% assigning the gains
bandGain(swir2Region) = scaleFactorMean(7);
bandGain(swir1Region) = scaleFactorMean(6);
bandGain(nirRegion) = scaleFactorMean(5);
bandGain(RedRegion) = scaleFactorMean(4);
bandGain(CARegion) = scaleFactorMean(2);

% selecting the starting point and ending  point
startpoint = blueGreenRegion(1);
endpoint = blueGreenRegion(end);
xValueInterpolationBlueGreenRegion = [hyperionWavelengths(startpoint - 1), modisCenterWavelength(2), ...
    hyperionWavelengths(endpoint+1)];
yValueInterpolationBlueGreenRegion = [bandGain(startpoint - 1),scaleFactorMean(3),bandGain(endpoint+1)];

interpolateWavelength = hyperionWavelengths(blueGreenRegion);
interpolateWavelengthGain = interp1(xValueInterpolationBlueGreenRegion, ...
    yValueInterpolationBlueGreenRegion,interpolateWavelength,'linear');

% inserting the intepolated data to bandGain Variabler
bandGain(blueGreenRegion) = interpolateWavelengthGain;
gainFactor = bandGain;    

% applying the gainFactor to the average of coincident pairs
gainCorrectedAverageHyperionReflectance = gainFactor.*averageHyperionReflectanceCoincident;

%plotting the average gain ccorrected Reflectance 
figure
plot(hyperionWavelengths,gainCorrectedAverageHyperionReflectance);
grid on
title('Average Hyperion Spectrum and banded Gain')
xlabel('Reflectance')
ylabel('wavelength(nm)');
hold on
plot(modisCenterWavelength,mean(bandedHyperionReflectance(:,2:end)),'lineStyle','None','Marker',...
    '^');

figure
plot(modisCenterWavelength,scaleFactorMean(2:end),'lineStyle','None','Marker',...
    '^');
title('Banded Gain Vs Wavelength')
xlabel('Banded Gain Factor')
ylabel('Wavelength')


figure
plot(hyperionWavelengths,bandGain,'LineStyle','None','MarKer','.','MarkerSize',30);
title('Gain Factor Vs Wavelength')
xlabel('Gain Factor')
ylabel('Wavelength')
ax = gca;
ax.FontSize = 20;
ax.YLim = [0.6  1.2];
hold on
plot(modisCenterWavelength,scaleFactorMean(2:end),'lineStyle','None','Marker',...
    '^');
grid on

% clearing the variables from the memory
% clearvars -except summaryAcquisitionDate coincidentPairsAngleSummary gainFactor ...
%      hyperionWavelengths scaleFactorDataSet averageHyperionReflectanceCoincident ...
%      gainCorrectedAverageHyperionReflectance
    
 save('Libya4GainCorrectionTest23.mat', 'gainFactor', 'coincidentPairsAngleSummary' ...
     ,'scaleFactorDataSet','summaryAcquisitionDate','hyperionWavelengths', ... 
     'averageHyperionReflectanceCoincident');
