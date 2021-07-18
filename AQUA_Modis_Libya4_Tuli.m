% AQUA MODIS collection 6.1 temporal trend analysis using data from NASA MCST Team (Amit Angal);
% Libya 4 temporal TOA Reflectance trend analysis;
clc;
clear all;
close all;
%Read Data files for Libya 4 from Excel file;
filename='Libya4_ONLY_AQUA_MODIS.xlsx';
AQUA_Modis_All_data_Libya4 = xlsread(filename);

% Extract TOA Reflectance information from column 8 to 13 (blue, green, red, NIR, SWIR 1 and SWIR 2 respectively)
AQUA_ModisLibya4.TOA_Reflectance = AQUA_Modis_All_data_Libya4 (:,8:13);

% as band 6 of AQUA has problem due to inoperable detectors near that band;
AQUA_ModisLibya4.TOA_Reflectance (:,5) = AQUA_ModisLibya4.TOA_Reflectance (:,6);
AQUA_ModisLibya4.TOA_Reflectance (:,6) = [];

% Year information of data aquisition
AQUA_ModisLibya4.Year = AQUA_Modis_All_data_Libya4 (:,2);

% Extract 4 angle information (order is SZA, SAA, VZA, VAA )
AQUA_ModisLibya4.Angle = AQUA_Modis_All_data_Libya4 (:,20:23);

% Extract standard deviation of each TOA reflectance from mean TOA reflectance (blue, green, red, NIR, SWIR 1 and SWIR 2 respectively)
AQUA_ModisLibya4.STD_MeanTOA = AQUA_Modis_All_data_Libya4 (:,14:19);

% SWIR 1 band is not appropriate, discard it from further works
AQUA_ModisLibya4.STD_MeanTOA (:,5) = AQUA_ModisLibya4.STD_MeanTOA (:,6);
AQUA_ModisLibya4.STD_MeanTOA (:,6) = [];


% Extract DSL and DOY information
AQUA_ModisLibya4.DSL =  AQUA_Modis_All_data_Libya4 (:,1);
AQUA_ModisLibya4.DOY = AQUA_Modis_All_data_Libya4 (:,3);

% Decimal Year Calculation from Year and DOY information
AQUA_ModisLibya4.Decimal_Year = AQUA_ModisLibya4.Year + AQUA_ModisLibya4.DOY/365;

% Modis_Data = DateInfoExtractDat('D:\sudan 1_Modis_Terra\Sudan1_ONLY_MODIS_forSDSU.dat');

% Suspected Outliers (cloudy scenes) find by 2 sigma process; It can keep
% information of outlier index for only last iteration-band 6
for b = 1:5
    outlier = abs(AQUA_ModisLibya4.TOA_Reflectance(:,b) - mean(AQUA_ModisLibya4.TOA_Reflectance(:,b))) > 2*std(AQUA_ModisLibya4.TOA_Reflectance(:,b));
    outlier_index = find(outlier);
end

cloudyScenes = outlier_index ; % they are really wierd, reflectance so high ; used band 6 (column 5 - SWIR 2) outlier indices 

% Exclude the Cloudy scenes from all the infomations
AQUA_ModisLibya4.TOA_Reflectance(cloudyScenes,:) = []; AQUA_ModisLibya4.STD_MeanTOA(cloudyScenes,:) = [];
AQUA_ModisLibya4.DOY(cloudyScenes,:) = []; AQUA_ModisLibya4.Angle(cloudyScenes,:) = []; AQUA_ModisLibya4.Year(cloudyScenes,:) = [];
AQUA_ModisLibya4.DSL(cloudyScenes,:) = [];AQUA_ModisLibya4.Decimal_Year(cloudyScenes,:) = [];

% Spatial Uncertainty calculation for making a threshold for outlier rejection
AQUA_ModisLibya4.spatial_uncertainty = (AQUA_ModisLibya4.STD_MeanTOA./AQUA_ModisLibya4.TOA_Reflectance)*100;
outlierIndex_wrt_spatial_U = find (AQUA_Modi sLibya4.spatial_uncertainty(:,5) >= 8 );

%{ when a scene has problem, it should show this for all bands less or more, that's why threshold was made w.r.t. band 1}

% Exclude the scenes which have spatial uncertainty greater than the threshold of 10 %
AQUA_ModisLibya4.TOA_Reflectance(outlierIndex_wrt_spatial_U,:) = []; AQUA_ModisLibya4.STD_MeanTOA(outlierIndex_wrt_spatial_U,:) = [];
AQUA_ModisLibya4.DOY(outlierIndex_wrt_spatial_U,:) = []; AQUA_ModisLibya4.Angle(outlierIndex_wrt_spatial_U,:) = []; 
AQUA_ModisLibya4.Year(outlierIndex_wrt_spatial_U,:) = [];AQUA_ModisLibya4.DSL(outlierIndex_wrt_spatial_U,:) = [];
AQUA_ModisLibya4.Decimal_Year(outlierIndex_wrt_spatial_U,:) = [];AQUA_ModisLibya4.spatial_uncertainty(outlierIndex_wrt_spatial_U,:) = [];

% Excluding scenes as cloudy from the website
cloudy = [1028,1015,952,736,733,731,675,664,551,520,462,461,389,332,329,268,198,184];

AQUA_ModisLibya4.TOA_Reflectance(cloudy,:) = []; AQUA_ModisLibya4.STD_MeanTOA(cloudy,:) = [];
AQUA_ModisLibya4.DOY(cloudy,:) = []; AQUA_ModisLibya4.Angle(cloudy,:) = []; 
AQUA_ModisLibya4.Year(cloudy,:) = [];AQUA_ModisLibya4.DSL(cloudy,:) = [];
AQUA_ModisLibya4.Decimal_Year(cloudy,:) = [];AQUA_ModisLibya4.spatial_uncertainty(cloudy,:) = [];


% Calculation of Temporal Uncertainty from Temporal Mean and Standard deviation information
for band = 1:5
    temporal_meanTOA_Reflectance(band) = mean (AQUA_ModisLibya4.TOA_Reflectance(:,band));
    temporal_STD_TOA_Reflectance (band)= std (AQUA_ModisLibya4.TOA_Reflectance(:,band));
    temporal_uncertainty_TOA_Reflectance (band) = (temporal_STD_TOA_Reflectance/temporal_meanTOA_Reflectance)*100;  
end

% 4 angle Linear BRDF normalization for reducing seasonality from the TOA reflectance trend

% keep angle informations in the corresponding variables
AQUA_ModisLibya4.SZA= AQUA_ModisLibya4.Angle(:,1);
AQUA_ModisLibya4.SAA = AQUA_ModisLibya4.Angle(:,2);
AQUA_ModisLibya4.VZA= AQUA_ModisLibya4.Angle(:,3);
AQUA_ModisLibya4.VAA = AQUA_ModisLibya4.Angle(:,4);

% Independent variable create by Spherical coordinates to rectangular coordinates conversion of the Solar and Azimuth angles
y1 = sind(AQUA_ModisLibya4.SZA).*sind(AQUA_ModisLibya4.SAA);
x1 = sind(AQUA_ModisLibya4.SZA).*cosd(AQUA_ModisLibya4.SAA);
y2 = sind(AQUA_ModisLibya4.VZA).*sind(AQUA_ModisLibya4.VAA);
x2 = sind(AQUA_ModisLibya4.VZA).*cosd(AQUA_ModisLibya4.VAA);

reflectance = AQUA_ModisLibya4.TOA_Reflectance;
n = length(reflectance);

for i = 1:5
    tbl = table(reflectance(:,i),y1,x1,y2,x2,'VariableNames',{'R' 'u' 'v' 'w' 'x'});
    lm = fitlm(tbl,'R~u+v+w+x');
    % sse(i)=lm.SSE;
    b = lm.Coefficients.Estimate;
    
    % same reference angle information used for all sensors BRDF correction
    y11 = sind(33).*sind(135); % reference chosen equal or close to mean
    x11 = sind(33).*cosd(135);
    y12 = sind(13).*sind(10);
    x12 = sind(13).*cosd(10);
    
    % ref_reflc(i)= b(1)+b(2)*y11+b(3)*x11+b(4)*y12+b(5)*x12; % reference reflectance value;
    ref_reflc(i) = mean (reflectance (:,i));
    
    for j = 1:n
        predicted_reflc(j,i)= b(1)+b(2)*y1(j)+b(3)*x1(j)+b(4)*y2(j)+b(5)*x2(j); % predicted reflectnace value:
        AQUA_ModisLibya4.corrected_reflc(j,i)= reflectance(j,i)*ref_reflc(i)/predicted_reflc(j,i); % corrected reflectance value;
        
    end
end
% BRDF uncertainty calculation (percentage)
AQUA_ModisLibya4.Aqua_brdfUncertainty = abs(((reflectance-predicted_reflc)./reflectance)*100);

% Temporal Uncertainty calculation of the mean TOA reflectance after BRDF
for band = 1:5
    temporal_uncertainty_after_BRDF (band) = (std(AQUA_ModisLibya4.corrected_reflc(:,band))/mean(AQUA_ModisLibya4.corrected_reflc(:,band)))*100;    
end

% plot the TOA reflectance before and after BRDF normalization individually band wise;

COL = {'c','g','r','m','k'};
col = {'[0.5 0.8 0.2]', '[0.2 0.6 0.8]', '[0.9 0.2 0.6]','[0.3 0.6 0.1]','[0.8 0.8 0.2]'};
Bands = {'Blue','Green','Red','NIR','SWIR2'};

for band = 1:5
    figure(band)
    plot(AQUA_ModisLibya4.Decimal_Year,AQUA_ModisLibya4.TOA_Reflectance(:,band),'LineStyle','none','Color',COL{band},'MarkerSize',20,'Marker','h',...
        'MarkerFaceColor',COL{band},'MarkerEdgeColor','black')
    xlabel('Decimal Year','FontSize',16,'FontWeight','bold','Color','k')
    ylabel('TOA Reflectance','FontSize',16,'FontWeight','bold','Color','k')
    ylim ([mean(AQUA_ModisLibya4.TOA_Reflectance(:,band)-0.25)  mean(AQUA_ModisLibya4.TOA_Reflectance(:,band)+0.25)])
    xlim ([2002 2019])
    title(['TOA Reflectance trend before and after BRDF normalization for ' Bands{band} '-' 'Band'],'FontSize',20);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 30); hold on; grid on; set (gcf,'color','white');
    plot(AQUA_ModisLibya4.Decimal_Year,AQUA_ModisLibya4.corrected_reflc(:,band),'LineStyle','none','Color',col{band},'MarkerSize',20,'Marker','o',...
        'MarkerFaceColor',col{band},'MarkerEdgeColor','black')
    legend({'Before Correction' 'After Correction' })  
end

%plot of the BRDF corrected meanTOA reflectance with spatial uncertainty bar for all 5 bands together; except band 5 - SWIR 1 band
for band = 1:5
    figure(band+5)
    errorbar(AQUA_ModisLibya4.Decimal_Year,AQUA_ModisLibya4.corrected_reflc(:,band),AQUA_ModisLibya4.spatial_uncertainty (:,band),...
        'o', 'MarkerSize', 16, 'MarkerFaceColor',COL{band},'MarkerEdgeColor','black','LineWidth',1.5,'Color',COL{band});
    xlabel('Decimal Year','FontSize',16,'FontWeight','bold','Color','k')
    ylabel('TOA Reflectance','FontSize',16,'FontWeight','bold','Color','k')
    xlim([2002 2019])
    title(['TOA Reflectance trend with spatial uncertainty for' '-' Bands{band}] ,'FontSize',16);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 30, 'Fontweight','bold');
    grid on; set (gcf,'color','white');
end

% Shapiro-Wilk test of normality of the lifetime meanTOA reflectance after BRDF using swtest function from file exchange
for band = 1:5
    [H(:,band), pValue(:,band), W(:,band)] =  swtest(AQUA_ModisLibya4.corrected_reflc(:,band),0.05); % Shapiro wilk test
    h(:,band) = adtest(AQUA_ModisLibya4.corrected_reflc(:,band)); %%%%%% Anderson-Darling test
    %%% Linear relationship check between decimal year and TOA reflectance
    correlation_matrix = corrcoef(AQUA_ModisLibya4.Decimal_Year(:,1),AQUA_ModisLibya4.corrected_reflc(:,band));
    % plot of normality check by histogram
    figure (band+10)
    histogram (AQUA_ModisLibya4.corrected_reflc(:,band),20,'Edgecolor','black','LineWidth',1)
    xlabel('Value of TOA Reflectance', 'FontSize',20,'FontWeight','bold','Color','k');
    ylabel('Frequency of Occurrence for TOA Reflectance','FontSize',20,'FontWeight','bold','Color','k');
    title(['Histogram of TOA Reflectance' '-' Bands{band} '-' 'Band'],'FontSize',20,'FontWeight','bold','Color','k');
    set(gca, 'FontSize', 20)
end

%% Near Coincident pair estimation of AQUA Modis with L8 - OLI
% process start time
startTime = datetime('now');

% defining the difference of 1 days between the images
acquistionTimeDifference = days((1));

% Libya 4 dat file doesn't have valid data to line 24 so skipping the data till 24
skipLine = 24;
fileNameAquaModis  = 'D:\AQUA_Libya 4\Libya4_ONLY_AQUA_MODIS.dat';
directoryOfL8 = 'Z:\ImageDrive\OLI.TIRS\L8\P181\R040\';
L8Product = 'LC1';

% reading the mean reflectance, solar and sensor angle
[acquistionDateAquaModis] = DateInfoExtractDat(fileNameAquaModis);

% reading the date information of the L8 Data
[acquistionDateL8] = DateInfoExtractMTL_L8(directoryOfL8,L8Product);

% sizes of the row
[L8yrow,~] = size(acquistionDateL8);
[Mrow,~] = size(acquistionDateAquaModis);

% cummulative addition of the DOY to get the DSL
dsl = dataset({cumsum(acquistionDateL8.DOY) 'DSL'});

% concatenating the dataset
acquistionDateL8 = [acquistionDateL8 dsl];

% OLI images will always be less than the Modis. selecting the L8 OLI time acquistion
L8DateSel = transpose(acquistionDateL8.Time);

% repeating the L8-OLI images
L8repeatDate = repmat((L8DateSel),Mrow,1);

%repeating the Modis images
AquaModisrepeatDate  = repmat(acquistionDateAquaModis.Time,1,L8yrow);

%taking the Date subraction to find the near coincident pairs
timeDifference  = AquaModisrepeatDate - L8repeatDate;

% elspsed duration in units of days by changing the format
timeDifference.Format = 'd';

% taking the absolute of all the time
timeDifference = abs(timeDifference);

% finding the days less than 1 days (overpass time difference is 2 hours and 39 minutes around)
[posRowCo,posColCo] = find(timeDifference < acquistionTimeDifference);

% finding the date of acqusition using unique to remove the repeated dates in the matrix dates on the L8
for i = 1:size(posRowCo)
    coincidentDatesL8(i) = (L8repeatDate(posRowCo(i),posColCo(i)));
    
    % finding the  match with the Hyperion
    match = find(coincidentDatesL8(i)== acquistionDateL8.Time);
    
    % finding the Date string
    L8DateString(i) = acquistionDateL8.DATE_STRING(match);
    
    % finding the DSL
    L8DSL(i) = acquistionDateL8.DSL(match);
    
    %finding the match
    L8RowCoSel(i) = match;
    
end
% dates on the Modis
for i = 1:size(posRowCo)
    coincidentDatesAquaModis(i) = (AquaModisrepeatDate(posRowCo(i),posColCo(i)));
    
    % finding the  match with the Modis
    match = find(coincidentDatesAquaModis(i) == acquistionDateAquaModis.Time);
    
    % finding the DSL
    modisDateString(i) = acquistionDateAquaModis.DATE_STRING(match);
    %finding the match
    modisRowCoSel(i) = match;
end
% date Difference between the dates obtained
differenceBetween  = abs(coincidentDatesL8- coincidentDatesAquaModis);

% process ending time
endTime = datetime('now');

% putting all the information in the table as summary table
summaryAcquisitionDate = table(coincidentDatesL8',coincidentDatesAquaModis',...
    differenceBetween',(L8DateString'),(modisDateString'),L8RowCoSel',...
    modisRowCoSel',(L8DSL)');

summaryAcquisitionDate.Properties.VariableNames = {'L8_Acquisition_Libya_4',...
    'Aqua_Modis_Acquisition_Libya_4','Time_Difference','L8_Date_String',...
    'Aqua_Modis_Date_String','L8_Row','Aqua_Modis_Row','L8DSL'};


%% Band Adjustment factor determination among the near coincident pairs of Aqua Modis and L8 OLI

% import the decimal years for those scenes (L8 and Modis for this script)
load ('D:\AQUA_Libya 4\OutlierFree_L8_Libya4_data_Alll');
%load ('D:\AQUA_Libya 4\Aqua_Modis_outlierFree_all_data_Libya4');

% Finding out decimal year from the data table of near coincident pairs of L8 and Aqua Modis using the similar formula as I used for my processed scene's
% decimal year estimation

dateFormat = 'yyyy-mm-dd';
DateString_L8 = datestr(summaryAcquisitionDate.L8_Acquisition_Libya_4,dateFormat);
DateString_Aqua_Modis = datestr(summaryAcquisitionDate.Aqua_Modis_Acquisition_Libya_4,dateFormat);
[doyL8,fractionL8] = date2doy(datenum(DateString_L8));
[doyMod,fractionMod] = date2doy(datenum(DateString_Aqua_Modis));

DateVector_L8 = datevec(DateString_L8,dateFormat);
Year_L8 = DateVector_L8(:,1);
near_coincidentPair_decimalYear_L8 = Year_L8 + doyL8/365;

DateVector_Modis = datevec(DateString_Aqua_Modis,dateFormat);
Year_AquaModis = DateVector_L8(:,1);
near_coincidentPair_decimalYear_AquaModis = Year_AquaModis + doyMod/365;


near_conincident_idx_L8 = ismember (Decimal_Year_L8,near_coincidentPair_decimalYear_L8);
rows_near_coincident_L8 = find (near_conincident_idx_L8);

near_conincident_idx_AquaModis = ismember (AQUA_ModisLibya4.Decimal_Year, near_coincidentPair_decimalYear_AquaModis);
rows_near_coincident_AquaModis = find (near_conincident_idx_AquaModis);

% matched near coincident pair indices from Aqua Modis
final_rows_near_coincident_AquaModis = rows_near_coincident_AquaModis(1:49,:);

% extracting the processed reflectance from Aqua Modis only for the near conincident indices
near_coincident_reflcetance_AquaModis = AQUA_ModisLibya4.corrected_reflc (final_rows_near_coincident_AquaModis,:);

% extracting the processed reflectance from L8 only for the near coincident indices (BRDF corrected meanTOA reflectance)
corrected_reflc (:,1) = [];
corrected_reflc (:,5) = [];
near_coincident_reflectance_L8 = corrected_reflc(rows_near_coincident_L8,:);

% ratio of the reflectance from L8 OLI and Aqua Modis for near conincident pairs only
ratioReflectance_nearcoincidentPairs = near_coincident_reflectance_L8./near_coincident_reflcetance_AquaModis;

% Band adjustment factor is the average of the ratio of reflectance of near conincident pairs
bAdj_factor = mean (ratioReflectance_nearcoincidentPairs);

% to reduce reflectance difference after all type of normalization though the use of band adjustment factor
bandAdjusted_TOA_AquaModis = AQUA_ModisLibya4.corrected_reflc.*bAdj_factor;

% band adjustment uncertainty calculation for Aqua MODIS (single uncertainty for all the TOA reflectance);
standardDeviation_bAdj_AquaModis = std(ratioReflectance_nearcoincidentPairs);
Aqua_bandAdjUncertainty = (standardDeviation_bAdj_AquaModis./bAdj_factor)*100;

% percentage difference of mean of meanTOA L7 before and after Band ajustment from L8 OLI
percent_diff_before_BAdj_AquaModis = (mean(AQUA_ModisLibya4.corrected_reflc) - mean(corrected_reflc))./ mean(corrected_reflc)*100
percent_diff_after_BAdj_AquaModis  = (mean(bandAdjusted_TOA_AquaModis) - mean(corrected_reflc))./ mean(corrected_reflc)*100


% specifying color for 6 different bands of two instruments
col = {'c','g','r','m','k'};
COL = {'[0.5 0.8 0.2]', '[0.2 0.6 0.8]', '[0.9 0.2 0.6]','[0.3 0.6 0.1]','[0.8 0.8 0.2]'};

% plot of L8 OLI and Aqua MOdis TOA trend after BRDF without band adjustment process
for b = 1:5
    figure (b+15)
    plot (AQUA_ModisLibya4.Decimal_Year,AQUA_ModisLibya4.corrected_reflc (:,b),'LineStyle','none','Color',....
        col{b},'MarkerSize',20,'Marker','o','MarkerFaceColor',col{b},'MarkerEdgeColor','black','LineWidth',1) % Aqua Modis reflectance
    ylim ([mean(AQUA_ModisLibya4.corrected_reflc(:,b)-0.20)   mean(AQUA_ModisLibya4.corrected_reflc(:,b)+0.20)]);
    xlim ([2002 2019])
    xlabel('Decimal Year','FontSize',16,'FontWeight','bold','Color','k')
    ylabel('TOA Reflectance','FontSize',16,'FontWeight','bold','Color','k')
    title(['BRDF corrected TOA Reflectance trend for ' Bands{b} '-' 'Band'],'FontSize',12);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 30)
    grid on; hold on;set (gcf,'color','white');
    plot (Decimal_Year_L8,corrected_reflc (:,b),'LineStyle','none','Color',....
        COL{b},'MarkerSize',20,'Marker','o','MarkerFaceColor',COL{b}, 'MarkerEdgeColor','black','LineWidth',1) % L8 OLI reflectance
    legend ('Aqua Modis', 'L8 OLI')
end



% plot of L8 OLI and Aqua MOdis TOA trend after BRDF and Band adjustment
for b = 1:5
    figure (b+20)
    plot (AQUA_ModisLibya4.Decimal_Year,bandAdjusted_TOA_AquaModis (:,b),'LineStyle','none','Color',....
        col{b},'MarkerSize',20,'Marker','o','MarkerFaceColor',col{b},'MarkerEdgeColor','black','lineWidth',1) % Aqua Modis reflectance
    ylim ([mean(bandAdjusted_TOA_AquaModis (:,b)-0.20)   mean(bandAdjusted_TOA_AquaModis (:,b)+0.20)]);
    xlim ([2002 2019])
    xlabel('Decimal Year','FontSize',16,'FontWeight','bold','Color','k')
    ylabel('TOA Reflectance','FontSize',16,'FontWeight','bold','Color','k')
    title(['BRDF corrected and band adjusted TOA Reflectance trend for ' Bands{b} '-' 'Band'],'FontSize',12);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 30)
    grid on; hold on;set (gcf,'color','white');
    plot (Decimal_Year_L8,corrected_reflc (:,b),'LineStyle','none','Color',....
        COL{b},'MarkerSize',20,'Marker','o','MarkerFaceColor',COL{b}, 'MarkerEdgeColor','black','lineWidth',1) % L8 OLI reflectance
    legend ('Aqua Modis', 'L8 OLI')
end

% --------------------------------------Harmonization--------------------------------------------------------------------------
% combining TOA data from all 5 sensors including new one AQUA Modis

% import data processed for other 4 sensors rather than Aqua in different folder
load ('Outlier_Free_All_data_Libya4_L7_ETM+');
load ('Outlier_free_data_Modis_Libya4');
load ('OutlierFree_L8_Libya4_data_Alll');
load('OutlierFree_All_Data_S2A_Libya4');

%import band adjustment and uncertainty information for these 4 instruments
load ('L7_L8_S2A_Terra_BandAdjInfo_All4typesUncertainty');

% concatenate all the Decimal years from 5 instruments
d1 = Decimal_Year_L7;
d2 = Decimal_Year_Modis; % Modis-Terra
d3 = AQUA_ModisLibya4.Decimal_Year;
d4 = Decimal_Year_L8;
d5 = Decimal_Year_S2A;
decimalYear_all_sensors = vertcat (d1,d2,d3,d4,d5);

% cocatenate all individual BRDF normalized band adjusted TOA reflectance from 5 different sensors to do statistical test
corrected_reflc_L8 = corrected_reflc(:,2:7);

corrected_reflc_L8 (:,5) = [];                   % remove SWIR 1 band;
bandAdj_TOA_L7 (:,5) = [];                       % remove SWIR 1 band;
bandAdjusted_TOA_Modis (:,5) = [];               % remove SWIR 1 band;
bandAdj_TOA_S2A (:,5) = [];                      % remove SWIR 1 band;

BRDF_bAdjusted_TOA_All = vertcat (bandAdj_TOA_L7,bandAdjusted_TOA_Modis,bandAdjusted_TOA_AquaModis,corrected_reflc_L8,bandAdj_TOA_S2A);

COL={'c','g','r','m','k'};
Bands = {'Blue', 'Green', 'Red', 'NIR','SWIR 2'};

for b = 1:1:5
    figure (100)
    plot(decimalYear_all_sensors,BRDF_bAdjusted_TOA_All(:,b),'LineStyle','none','Color',....
        COL{b},'MarkerSize',20,'Marker','o','MarkerFaceColor',COL{b},'MarkerEdgeColor','k','lineWidth',0.7)
    ylim ([0.15 0.75])
    xlim ([1999 2019])
    xlabel('Decimal Year','FontSize',16,'FontWeight','bold','Color','k')
    ylabel('TOA Reflectance','FontSize',16,'FontWeight','bold','Color','k')
    title(['5 Sensors harmonized Combined TOA Reflectance trend of' ' ' 'Libya 4'],'FontSize',20);
    xt = get(gca, 'XTick');set(gca, 'FontSize', 30, 'fontWeight', 'bold'); set (gcf, 'color','white'); grid on;hold on;
    % least square regression fit in data
    fit = polyfit(decimalYear_all_sensors,BRDF_bAdjusted_TOA_All(:,b),1);
    yfit = polyval(fit,decimalYear_all_sensors); hold on
    plot(decimalYear_all_sensors,yfit,'b-','MarkerSize',15,'LineWidth',0.7);hold on
    legend ({'Blue','Fitted line Blue','Green','Fitted line Green','Red','Fitted line Red','NIR','Fitted line NIR',...
        'SWIR 2','Fitted line SWIR 2'}, 'Location','northwestoutside');
    
end

% combine dsl(day since launch information) from 5 instruments
dsl_L7 = dsl_vect_L7; dsl_Modis = DSL_Modis_Libya4; dsl_Aqua = AQUA_ModisLibya4.DSL; dsl_L8 = dsl_vect_L8_Libya4; dsl_S2A = dsl_vect_S2A;
dsl_all_sensors = vertcat (dsl_L7,dsl_Modis,dsl_Aqua,dsl_L8,dsl_S2A);

% Check the distribution of the combined data and check if there is any linear relationship between TOA reflectance and decimal year
% Shapiro-Wilk test of normality of the lifetime combined meanTOA reflectance after BRDF and band adjustment using swtest function from file exchange

bandNames = {'Blue','Green','Red','NIR','SWIR2'};
for band = 1:5
    [H(:,band), pValue(:,band), W(:,band)] =  swtest(BRDF_bAdjusted_TOA_All(:,band),0.05); % Shapiro wilk test
    h(:,band) = adtest(BRDF_bAdjusted_TOA_All(:,band)); %%%%%% Anderson-Darling test
    % Linear relationship check between decimal year and TOA reflectance;
    % if dsl is also used to check linear relationship, for some bands like blue,green and swir 2 linearity failed
    [correlation_matrix,P_Linearity] = corrcoef(dsl_all_sensors(:,1),BRDF_bAdjusted_TOA_All(:,band))
    % plot of normality check by histogram
    figure(band+13)
    histogram (BRDF_bAdjusted_TOA_All(:,band),20,'Edgecolor','black','LineWidth',1.5)
    xlabel('Value of TOA Reflectance', 'FontSize',20,'FontWeight','bold');
    ylabel('Frequency of Occurrence for TOA Reflectance','FontSize',20,'FontWeight','bold');
    title(['Histogram of Harmonized TOA Reflectance data over Libya 4' '-' bandNames{band} '-' 'Band'],'FontSize',25);
    set(gca, 'FontSize', 25, 'FontWeight','bold');set (gcf, 'color','white'); grid on;
    % check autocorrelation
    figure (band+18)
    autocorr(BRDF_bAdjusted_TOA_All(:,band),50,[],3) % show autocorrelation function for 1st 50 lags
    % conduct the Ljung-Box Q test for autocorrelation because need this information in Seasonal mann kendall test p value slection
    [h_autocorrelation,p,Qstat,crit] = lbqtest(BRDF_bAdjusted_TOA_All(:,band),'Lags',[12,24,36,48])
end


%------------- T test with respect to decimal years/dsl---------------------
for b = 1:5
    T = table(decimalYear_all_sensors,BRDF_bAdjusted_TOA_All(:,b),'VariableNames',{'Decimal_Year' 'Reflectance' });
    mdl1 = fitlm(T);
    slope_decimalYears(:,b)= mdl1.Coefficients.Estimate(2);
    SE_Slope_decimalYears(:,b)= mdl1.Coefficients.SE(2);
    P_value_decimalYears(:,b)= mdl1.Coefficients.pValue(2);
    t_value_decimalYears(:,b)= mdl1.Coefficients.tStat(2);   
end

for b = 1:5
    T = table(dsl_all_sensors,BRDF_bAdjusted_TOA_All(:,b),'VariableNames',{'Decimal_Year' 'Reflectance' });
    mdl2 = fitlm(T);
    slope_dsl(:,b)= mdl2.Coefficients.Estimate(2);
    SE_Slope_dsl(:,b)= mdl2.Coefficients.SE(2);
    P_value_dsl(:,b)= mdl2.Coefficients.pValue(2);
    t_value_dsl(:,b)= mdl2.Coefficients.tStat(2);   
end

   
%% Chi Square test considering all types of calculation uncertainty (Sensor uncertainty, Spatial Uncertainty, BRDF Uncertainty, Band Adjustment Uncertainty)
%Sensor radiometric calibration Uncertainty arrays
L7_sensorUncertainty = zeros(232,5); 
Terra_sensorUncertainty = repmat(2,291,5);
Aqua_sensorUncertainty = repmat(2,1124,5);
L8_sensorUncertainty = repmat(2,61,5);
S2A_sensorUncertainty = repmat(5,58,5);

% concatenate sensor uncertainty arrays
U_sensor = vertcat (L7_sensorUncertainty,Terra_sensorUncertainty,Aqua_sensorUncertainty,L8_sensorUncertainty,S2A_sensorUncertainty);

% Spatial Uncertainty of the calculated mean TOA reflectances sensor wise
L7_spatialUncertainty  = (TOA_reflectance_L7_STD./TOA_reflectance_L7)*100; L7_spatialUncertainty (:,5) = [];
Terra_spatialUncertainty = (STD_MeanTOA./TOA_M_Libya4_Modis)*100;Terra_spatialUncertainty (:,5) = [];
Aqua_spatialUncertainty = AQUA_ModisLibya4.spatial_uncertainty;
L8_spatialUncertainty = (STDTOA(:,2:7)./meanTOA(:,2:7))*100;L8_spatialUncertainty(:,5)=[];
S2A_spatialUncertainty = (STD_Refl_S2A(2:7)./Refl_S2A(:,2:7))*100;S2A_spatialUncertainty (:,5)=[];

% combine spatial uncertainty for all 4 sensors
U_spatial = vertcat (L7_spatialUncertainty,Terra_spatialUncertainty,Aqua_spatialUncertainty,L8_spatialUncertainty,S2A_spatialUncertainty);

% combine BRDF model calculation uncertainty for all 4 sensors;SWIR 1 band is being removed
L7_brdfUncertainty (:,5)=[];
Terra_brdfUncertainty (:,5)=[];
L8_brdfUncertainty = L8_brdfUncertainty (:,2:7) ; 
L8_brdfUncertainty (:,5)=[];

S2A_brdfUncertainty = S2A_brdfUncertainty (:,2:7) ; 
S2A_brdfUncertainty (:,5)=[];

U_BRDF = vertcat (L7_brdfUncertainty,Terra_brdfUncertainty,AQUA_ModisLibya4.Aqua_brdfUncertainty,L8_brdfUncertainty ,S2A_brdfUncertainty);

% Combine band adjustment factor calculation uncertainty for all 4 sensors
% As landsat 8 OLI doesn't need band adjustment; create an array with zero values for L8 while combining
L8_bandAdjUncertainty = zeros(61,5);
L7_bandAdjUncertainty (:,5) = [];
Terra_bandAdjUncertainty (:,5) = [];
Aqua_bandAdjUncertainty = repmat (Aqua_bandAdjUncertainty,1124,1)
S2A_bandAdjUncertainty (:,5) = [];

U_bandAdjustmentFactor = vertcat(L7_bandAdjUncertainty,Terra_bandAdjUncertainty,Aqua_bandAdjUncertainty,L8_bandAdjUncertainty,S2A_bandAdjUncertainty);

% total absolute uncertainty calculation from different uncertainties
total_Uncertainty = sqrt(U_sensor.^2+ U_spatial.^2 + U_BRDF.^2 + U_bandAdjustmentFactor.^2);
% absolute value of the total uncertainty (standard deviation)
sigma = (total_Uncertainty.*BRDF_bAdjusted_TOA_All)./100;


% Plotting lifetime TOA reflectance with uncertainty
bandNames={'Blue','Green','Red','NIR','SWIR2'};
COL = {'c','g','r','m','k'};
for j = 1:5
    figure(j+31)
    errorbar(decimalYear_all_sensors,BRDF_bAdjusted_TOA_All(:,j), sigma(:,j), 'o', 'MarkerSize', 20, 'MarkerFaceColor',COL{j},...
        'MarkerEdgeColor','black','LineWidth',1.5, 'Color',COL{j},'LineWidth',1.5);
    ylim([mean(BRDF_bAdjusted_TOA_All(:,j))-0.25 mean(BRDF_bAdjusted_TOA_All(:,j))+0.25])
    xlim([1999 2019])
    xlabel('Decimal Year','FontSize',16,'FontWeight','bold','Color','k')
    ylabel('mean TOA Reflectance','FontSize',16,'FontWeight','bold','Color','k')
    title(['Combined TOA Reflectance for Libya 4 with total uncertainty for' '-' bandNames{j}] ,'FontSize',16);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 30, 'Fontweight','bold');
    grid on; set (gcf,'color','white');
end


% Chi square test considering all the uncertainties in a single TOA reflectance. This test will help to find out which fit (with slope or only the mean) is the good one
% w.r.t chi squared values

for b = 1:5
    T = table(dsl_all_sensors(:,1),BRDF_bAdjusted_TOA_All(:,b),'VariableNames',{'DSL' 'Reflectance' });
    mdl3 = fitlm(T);
    disp(mdl3);
    slope_all_sensors(:,b)= mdl3.Coefficients.Estimate(2);
    Intercept_all_sensors(:,b)= mdl3.Coefficients.Estimate(1);
    SE_of_slope_all_sensors(:,b)= mdl3.Coefficients.SE(2);
    Expected_TOA_Reflc(:,b) = mdl3.Fitted;
    Mean(b) = mean(BRDF_bAdjusted_TOA_All(:,b));
    
    % Finding chi square values for the fitted line with slope
    X_fit1 =(((BRDF_bAdjusted_TOA_All(:,b)- Expected_TOA_Reflc)./sigma(:,b)).^2);
    Chi_Square_fit1 = sum(X_fit1);%%%%%%%Chi-square test for 1st fit
    n = size(dsl_all_sensors(:,1),1); %%%%sample size
    df = n-1;%%%%degrees of freedom
    Reduced_Chi_Square_fit1 = Chi_Square_fit1./df;
    
    % Finding Chi square values for fitted line without slope, considering only the average of TOA reflectance;
    X_fit2 =(((BRDF_bAdjusted_TOA_All(:,b)- Mean)./sigma(:,b)).^2);
    Chi_Square_fit2 = sum(X_fit2);%%%%%%%Chi-square test for 2nd fit;
    Reduced_Chi_Square_fit2 = Chi_Square_fit2./df; 
    
    % AIC (Akaike Information criteria) for getting best model between this two
    N = length (BRDF_bAdjusted_TOA_All);
    AIC_fit1 = Chi_Square_fit1 + 2*2 + 2*2*(2+1)/(N-2-1);
    AIC_fit2 = Chi_Square_fit2 + 2*1 + 2*1*(1+1)/(N-1-1);
end

% ----------------------------Plot after BRDF and Band adjutment compensation individually------------------------------------------------------------
for b = 1:1:5
    figure (b+25)
    
    % Aqua-modis plot
    plot(AQUA_ModisLibya4.Decimal_Year,bandAdjusted_TOA_AquaModis(:,b),'LineStyle','none','Color',.....
        'green','MarkerSize',20,'Marker','o','MarkerFaceColor','m','MarkerEdgeColor','black','LineWidth',1.5);hold on
    ylim([mean(bandAdjusted_TOA_AquaModis(:,b))-0.20  mean(bandAdjusted_TOA_AquaModis(:,b))+0.20]);
    xlim ([1999 2019])
    xlabel('Decimal Year','FontSize',16,'FontWeight','bold','Color','k')
    ylabel('TOA Reflectance','FontSize',16,'FontWeight','bold','Color','k')
    title(['BRDF corrected band adjusted TOA Reflectance trend for ' Bands{b} '-' 'Band' ' ' '(Libya 4)'],'FontSize',12);
    xt = get(gca, 'XTick'); set(gca, 'FontSize', 30, 'FontWeight','bold'); set (gcf, 'color', 'white');grid on; hold on
    fitModis = polyfit (AQUA_ModisLibya4.Decimal_Year,bandAdjusted_TOA_AquaModis (:,b),1);
    yfitModis = polyval (fitModis,AQUA_ModisLibya4.Decimal_Year);
    hold on
    plot(AQUA_ModisLibya4.Decimal_Year,yfitModis,'m-','MarkerSize',15,'LineWidth',2.5);
    hold on
    
    % L7-ETM+ plot
    plot(Decimal_Year_L7,bandAdj_TOA_L7(:,b),'LineStyle','none','Color',....
        'blue','MarkerSize',20,'Marker','o','MarkerFaceColor','blue','MarkerEdgeColor','black','LineWidth',1.5)
    fitL7 = polyfit(Decimal_Year_L7,bandAdj_TOA_L7(:,b),1);
    yfit7 = polyval(fitL7,Decimal_Year_L7);
    hold on
    plot(Decimal_Year_L7,yfit7,'b-','MarkerSize',15,'LineWidth',2.5);
    hold on
    
    % MODIS-TERRA plot
    plot(Decimal_Year_Modis,bandAdjusted_TOA_Modis (:,b),'LineStyle','none','Color',.....
        'green','MarkerSize',20,'Marker','o','MarkerFaceColor','green','MarkerEdgeColor','black','LineWidth',1.5);hold on
    fitModis = polyfit (Decimal_Year_Modis, bandAdjusted_TOA_Modis(:,b),1);
    yfitModis = polyval (fitModis,Decimal_Year_Modis);
    hold on
    plot(Decimal_Year_Modis,yfitModis,'g-','MarkerSize',15,'LineWidth',2.5);
    hold on
       
    % Landsat 8-OLI plot;as L8 is reference for band adjustment it's raw filtered data used
    plot (Decimal_Year_L8,corrected_reflc_L8(:,b),'LineStyle','none','Color',....
        'red','MarkerSize',20,'Marker','o','MarkerFaceColor','red','MarkerEdgeColor','black','LineWidth',1.5); hold on
    fitL8 = polyfit (Decimal_Year_L8,corrected_reflc_L8(:,b),1);
    yiftL8 = polyval (fitL8,Decimal_Year_L8);
    hold on
    plot(Decimal_Year_L8,yiftL8,'r-','MarkerSize',15,'LineWidth',2.5);
    hold on
    
    % S2A-MSI plot
    plot (Decimal_Year_S2A,bandAdj_TOA_S2A (:,b),'LineStyle','none','Color',....
        'black','MarkerSize',20,'Marker','o','MarkerFaceColor','black','MarkerEdgeColor','yellow','LineWidth',1.5); hold on
    fitS2A = polyfit (Decimal_Year_S2A,bandAdj_TOA_S2A (:,b),1);
    yiftS2A = polyval (fitS2A,Decimal_Year_S2A);
    hold on
    plot(Decimal_Year_S2A,yiftS2A,'k-','MarkerSize',15,'LineWidth',2.5);
    legend ({'Aqua-Modis','Aqua-modis fitted line','L7-ETM+','L7-Fitted line', 'Terra-Modis','Modis-Fitted line', 'L8-OLI',...
        'L8-Fitted line', 'S2A-MSI', 'S2A-Fitted line', 'Location','northwest'}, 'FontSize',16, 'Orientation', 'vertical')
    
end






















