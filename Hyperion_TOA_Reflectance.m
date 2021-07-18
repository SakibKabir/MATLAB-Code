%clc
clear all
%close all

startTime = datetime('now');

if ispc 
%defining the directory 
    dirName = 'Z:\ImageDrive\Hyperion\EO1\P195\R045\';
elseif isunix
    dirName = '/home/sakib.kabir/zdrive/ImageDrive/Hyperion/EO1/P033/R037/';
end
% defining the product 
product = 'L1T';

%get the directory list
dir_list = dir(dirName);

% launch date of the hyperion
launchDate =  datetime('2000/11/21','InputFormat','yyyy/MM/dd');
launchdateNum = datenum(launchDate);

% removing the first two rows
dir_list = dir_list(~ismember({dir_list.name},{'.','..'}));

% number of images in the directory
numberOfImages = length(dir_list);

%equator= 10000000;
UL_x = 453169;
UL_y = 2142263;
UR_x = 459671;
UR_y = 2140674;
LR_x = 450560;
LR_y = 2100280;
LL_x = 443742;
LL_y = 2099793;

% % defining the ROI % most likly libya 4
% XMAP1 =  766818; ...759472;     %758630; %705400  %767930;
% YMAP1 =  3187123; ...3189887;     %3179320;  %3204520  %3200530;
% XMAP2 =  769718; ...775000;     %770350;  %793400;  %773320;
% YMAP2 =  3182806; ...3143865;    %3167420;% 3196631; %3194260;
% 
% % defininng new ROI
% XMAP1 = 778335;
% YMAP1 = 3237675;
% XMAP2 = 781335;
% YMAP2 = 3236475;
% 
% ROI.XMAP1 = 778335;
% ROI.XMAP2 = 781335;
% ROI.YMAP1 = 3237675;
% ROI.YMAP2 = 3236475;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% defining the ROI for Sudan 1 (Absolute calibration)
% XMAP1 = 614210;
% XMAP2 = 617835;
% YMAP1 = 2377968;
% YMAP2 = 2375569;
% 
% ROI.XMAP1 = 614210;
% ROI.XMAP2 = 617835;
% ROI.YMAP1 = 2377968;
% ROI.YMAP2 = 2375569;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%}

% %Niger 1 roi
% XMAP1 =  564100; ...759472;     %758630; %705400  %767930;
% YMAP1 =  2235245; ...3189887;     %3179320;  %3204520  %3200530;
% XMAP2 =  566909; ...775000;     %770350;  %793400;  %773320;
% YMAP2 =  2230899; ...3143865;    %3167420;% 3196631; %3194260;

% XMAP1 = 557200;
% YMAP1 = 2292220;
% XMAP2 = 581735;
% YMAP2 = 2212570;
% 
% ROI.XMAP1 = 557200;
% ROI.YMAP1 = 2292220;
% ROI.XMAP2 = 581735;
% ROI.YMAP2 = 2212570;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % defining the ROI for the Egypt 1 (Absolute Calibration Method)
% XMAP1 = 371190;
% XMAP2 = 375453;
% YMAP1 = 3073616;
% YMAP2 = 3065888;
% 
% ROI.XMAP1 = 371190;
% ROI.XMAP2 = 375453;
% ROI.YMAP1 = 3073616;
% ROI.YMAP2 = 3065888;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% defining the ROI for the Libya 1 (Absolute Calibration Method)
% XMAP1 = 334254;
% XMAP2 = 337978;
% YMAP1 = 2720530;
% YMAP2 = 2716380;
% 
% ROI.XMAP1 = 334254;
% ROI.XMAP2 = 337978;
% ROI.YMAP1 = 2720530;
% ROI.YMAP2 = 2716380;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % defining the ROI for the Niger 2 (Absolute Calibration Method)
% XMAP1 = 662392;
% XMAP2 = 666854;
% YMAP1 = 2367988;
% YMAP2 = 2360825;
% 
% ROI.XMAP1 = 662392;
% ROI.XMAP2 = 666854;
% ROI.YMAP1 = 2367988;
% ROI.YMAP2 = 2360825;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ESUN derivered from CHkur Solar Model
esun=[989.15
1100.38000000000
1096.10
1026.55
1241.03
1658.58
1724.74
1643.21
1740.79
1966.18
2037.54
2004.42
2015.33
1928.46
1933.40
1915.66
1830.66
1883.64
1875.99
1872.42
1843.78
1845.40
1843.22
1791.31
1763.11
1733.16
1703.57
1669.99
1631.41
1574.55
1537.62
1523.46
1472.65
1447.53
1397.29
1359.63
1263.93
1273.65
1283.29
1268.52
1231.76
1202.21
1190.14
1151.97
1127.86
1100.06
1068.93
1045.32
1015
968.970
962.650
958.140
951.040
944.040
944.600
933.570
929.640
919.370
913.520
893.950
888.68
875.97
869.83
847.66
843.20
838.31
832.34
815.11
811.16
799.62
793.69
784.08
779.84
770.64
765.96
756.79
751.88
743.08
737.90
723.36
717.78
711.22
705.81
697.91
693.43
683.01
677.01
668.55
663.65
655.47
648.75
637.88
621.11
605.72
586.54
580.94
572.40
560.23
549.10
540.26
532.37
519.04
510.64
499.26
490.21
479.65
476.44
469.19
460.72
452.05
443.90
435.50
425.30
412.66
413.12
405.60
397.04
391.62
385.10
379.45
371.59
365.03
359.33
354.47
348.22
341.17
334.41
326.58
322.80
316.72
312.20
306.98
300.73
295.47
291.77
284.70
284.49
280.31
275.98
271.51
266.30
260.80
253.18
247.41
247.71
242.71
238.67
237.89
228.45
226.44
224.92
218.70
210.18
209.62
207
201.61
197.65
190.64
185.03
185.22
182.80
179.56
175.56
172.90
169
164.26
157.81
158.66
155.30
152.51
148.99
142.07
140.06
139.41
137.95
136.53
133.81
130.24
124.87
124.01
123.97
121.58
119.25
117.62
115.81
114.31
111.63
109.26
107.85
106.12
103.71
102.31
100.43
98.340
97.210
95.480
93.600
92.370
90.880
89.230
85.170
84.900
84.370
83.440
81.730
80.560
79.370
78.100
76.690
75.320
74.130
73.100
71.640
70.180
69.430
68.240
66.460
65.770
65.040
63.280
62.810
61.63
60.11
59.86
59.11
57.45
57.03
56.24
55.10
54.09
53.72
52.68
51.70
51.42
50.65
49.83
49.36
48.42
47.61
46.94
46.35
45.81];

% hyperion wavelenght array
hyperionbandwavelengths=1.0e+003*[0.3556
    0.3658
    0.3759
    0.3861
    0.3963
    0.4065
    0.4166
    0.4268
    0.4370
    0.4472 
    0.4573 %11
    0.4675 %ETM+ band 1
    0.4777 %
    0.4879 %
    0.4980 %
    0.5082 %
    0.5184 %17
    0.5286 
    0.5387 %19
    0.5489 %
    0.5591 %
    0.5693 %ETM+ band 2
    0.5795 %
    0.5896 %
    0.5998 %
    0.6100 %26
    0.6201
    0.6303 %28
    0.6405 %
    0.6507 %
    0.6609 %ETM+ band 3
    0.6710 %
    0.6812 %33
    0.6914
    0.7016
    0.7117
    0.7219
    0.7321
    0.7422
    0.7524
    0.7626
    0.7728
    0.7830 %43
    0.7931 %
    0.8033 %
    0.8135 % 
    0.8236 %
    0.8338 %
    0.8440 % ETM+ band 4
    0.8542 %
    0.8644 %
    0.8745 %
    0.8847 %
    0.8949 %54
    0.9050
    0.9152
    0.9254
    0.9356
    0.9458
    0.9559
    0.9661
    0.9763
    0.9865
    0.9966
    1.0068
    1.0169
    1.0271
    1.0373
    1.0475
    1.0576
    0.8519
    0.8620
    0.8721
    0.8822
    0.8923
    0.9024
    0.9125
    0.9225
    0.9326
    0.9427
    0.9528
    0.9629
    0.9730
    0.9831
    0.9932
    1.0033
    1.0133
    1.0234
    1.0334
    1.0435
    1.0536
    1.0637
    1.0738
    1.0839
    1.0940
    1.1041
    1.1141
    1.1242
    1.1343
    1.1440
    1.1540
    1.1640
    1.1740
    1.1840
    1.1940
    1.2050
    1.2150
    1.2250
    1.2350
    1.2450
    1.2550
    1.2650
    1.2750
    1.2850
    1.2950
    1.3050
    1.3160
    1.3260
    1.3360
    1.3460
    1.3560
    1.3660
    1.3760
    1.3860
    1.3960
    1.4060
    1.4160
    1.4260
    1.4370
    1.4470
    1.4570
    1.4670
    1.4770
    1.4870
    1.4970
    1.5070
    1.5170
    1.5270
    1.5370
    1.5480
    1.5580 % 141
    1.5680 %
    1.5780 %
    1.5880 %
    1.5980 %
    1.6080 %
    1.6180 %
    1.6280 %
    1.6380 %
    1.6480 %
    1.6590 %
    1.6690 %%ETM+ band 5
    1.6790 %
    1.6890 %
    1.6990 %
    1.7090 %
    1.7190 % 
    1.7290 %
    1.7390 %
    1.7490 % 160
    1.7590
    1.7690
    1.7800
    1.7900
    1.8000
    1.8100
    1.8200
    1.8300
    1.8400
    1.8500
    1.8600
    1.8700
    1.8800
    1.8910
    1.9010
    1.9110
    1.9210
    1.9310
    1.9410
    1.9510
    1.9610
    1.9710
    1.9810
    1.9910
    2.0020
    2.0120
    2.0220
    2.0320
    2.0420
    2.0520
    2.0620
    2.0720
    2.0820
    2.0920 %194
    2.1020 %
    2.1130 %
    2.1230 %
    2.1330 %
    2.1430 %
    2.1530 %
    2.1630 %
    2.1730 %
    2.1830 %
    2.1930 %
    2.2030 %
    2.2130 %ETM+ band 7
    2.2240 %
    2.2340 %
    2.2440 %
    2.2540 %
    2.2640 %
    2.2740 %
    2.2840 %
    2.2940 %
    2.3040 %
    2.3140 %
    2.3240 %
    2.3350 %
    2.3450 %219
    2.3550
    2.3650
    2.3750
    2.3850
    2.3950
    2.4050
    2.4150
    2.4250
    2.4350
    2.4450
    2.4560
    2.4660
    2.4760
    2.4860
    2.4960
    2.5060
    2.5160
    2.5260
    2.5360
    2.5460
    2.5560
    2.5660
    2.5770];

% sundistance array
sundistance=[0.983310000000000;0.983300000000000;0.983300000000000;0.983300000000000;0.983300000000000;0.983320000000000;0.983330000000000;0.983350000000000;0.983380000000000;0.983410000000000;0.983450000000000;0.983490000000000;0.983540000000000;0.983590000000000;0.983650000000000;0.983710000000000;0.983780000000000;0.983850000000000;0.983930000000000;0.984010000000000;0.984100000000000;0.984190000000000;0.984280000000000;0.984390000000000;0.984490000000000;0.984600000000000;0.984720000000000;0.984840000000000;0.984960000000000;0.985090000000000;0.985230000000000;0.985360000000000;0.985510000000000;0.985650000000000;0.985800000000000;0.985960000000000;0.986120000000000;0.986280000000000;0.986450000000000;0.986620000000000;0.986800000000000;0.986980000000000;0.987170000000000;0.987350000000000;0.987550000000000;0.987740000000000;0.987940000000000;0.988140000000000;0.988350000000000;0.988560000000000;0.988770000000000;0.988990000000000;0.989210000000000;0.989440000000000;0.989660000000000;0.989890000000000;0.990120000000000;0.990360000000000;0.990600000000000;0.990840000000000;0.991080000000000;0.991330000000000;0.991580000000000;0.991830000000000;0.992080000000000;0.992340000000000;0.992600000000000;0.992860000000000;0.993120000000000;0.993390000000000;0.993650000000000;0.993920000000000;0.994190000000000;0.994460000000000;0.994740000000000;0.995010000000000;0.995290000000000;0.995560000000000;0.995840000000000;0.996120000000000;0.996400000000000;0.996690000000000;0.996970000000000;0.997250000000000;0.997540000000000;0.997820000000000;0.998110000000000;0.998400000000000;0.998680000000000;0.998970000000000;0.999260000000000;0.999540000000000;0.999830000000000;1.00012000000000;1.00041000000000;1.00069000000000;1.00098000000000;1.00127000000000;1.00155000000000;1.00184000000000;1.00212000000000;1.00240000000000;1.00269000000000;1.00297000000000;1.00325000000000;1.00353000000000;1.00381000000000;1.00409000000000;1.00437000000000;1.00464000000000;1.00492000000000;1.00519000000000;1.00546000000000;1.00573000000000;1.00600000000000;1.00626000000000;1.00653000000000;1.00679000000000;1.00705000000000;1.00731000000000;1.00756000000000;1.00781000000000;1.00806000000000;1.00831000000000;1.00856000000000;1.00880000000000;1.00904000000000;1.00928000000000;1.00952000000000;1.00975000000000;1.00998000000000;1.01020000000000;1.01043000000000;1.01065000000000;1.01087000000000;1.01108000000000;1.01129000000000;1.01150000000000;1.01170000000000;1.01191000000000;1.01210000000000;1.01230000000000;1.01249000000000;1.01267000000000;1.01286000000000;1.01304000000000;1.01321000000000;1.01338000000000;1.01355000000000;1.01371000000000;1.01387000000000;1.01403000000000;1.01418000000000;1.01433000000000;1.01447000000000;1.01461000000000;1.01475000000000;1.01488000000000;1.01500000000000;1.01513000000000;1.01524000000000;1.01536000000000;1.01547000000000;1.01557000000000;1.01567000000000;1.01577000000000;1.01586000000000;1.01595000000000;1.01603000000000;1.01610000000000;1.01618000000000;1.01625000000000;1.01631000000000;1.01637000000000;1.01642000000000;1.01647000000000;1.01652000000000;1.01656000000000;1.01659000000000;1.01662000000000;1.01665000000000;1.01667000000000;1.01668000000000;1.01670000000000;1.01670000000000;1.01670000000000;1.01670000000000;1.01669000000000;1.01668000000000;1.01666000000000;1.01664000000000;1.01661000000000;1.01658000000000;1.01655000000000;1.01650000000000;1.01646000000000;1.01641000000000;1.01635000000000;1.01629000000000;1.01623000000000;1.01616000000000;1.01609000000000;1.01601000000000;1.01592000000000;1.01584000000000;1.01575000000000;1.01565000000000;1.01555000000000;1.01544000000000;1.01533000000000;1.01522000000000;1.01510000000000;1.01497000000000;1.01485000000000;1.01471000000000;1.01458000000000;1.01444000000000;1.01429000000000;1.01414000000000;1.01399000000000;1.01383000000000;1.01367000000000;1.01351000000000;1.01334000000000;1.01317000000000;1.01299000000000;1.01281000000000;1.01263000000000;1.01244000000000;1.01223000000000;1.01205000000000;1.01186000000000;1.01165000000000;1.01145000000000;1.01124000000000;1.01103000000000;1.01081000000000;1.01060000000000;1.01037000000000;1.01015000000000;1.00992000000000;1.00969000000000;1.00946000000000;1.00922000000000;1.00898000000000;1.00874000000000;1.00850000000000;1.00825000000000;1.00800000000000;1.00775000000000;1.00750000000000;1.00724000000000;1.00698000000000;1.00672000000000;1.00646000000000;1.00620000000000;1.00593000000000;1.00566000000000;1.00539000000000;1.00512000000000;1.00485000000000;1.00457000000000;1.00430000000000;1.00402000000000;1.00374000000000;1.00346000000000;1.00318000000000;1.00290000000000;1.00262000000000;1.00234000000000;1.00205000000000;1.00177000000000;1.00148000000000;1.00119000000000;1.00091000000000;1.00062000000000;1.00033000000000;1.00005000000000;0.999760000000000;0.999470000000000;0.999180000000000;0.998900000000000;0.998610000000000;0.998320000000000;0.998040000000000;0.997750000000000;0.997470000000000;0.997180000000000;0.996900000000000;0.996620000000000;0.996340000000000;0.996050000000000;0.995770000000000;0.995500000000000;0.995220000000000;0.994940000000000;0.994670000000000;0.994400000000000;0.994120000000000;0.993850000000000;0.993590000000000;0.993320000000000;0.993060000000000;0.992790000000000;0.992530000000000;0.992280000000000;0.992020000000000;0.991770000000000;0.991520000000000;0.991270000000000;0.991020000000000;0.990780000000000;0.990540000000000;0.990300000000000;0.990070000000000;0.989830000000000;0.989610000000000;0.989380000000000;0.989160000000000;0.988940000000000;0.988720000000000;0.988510000000000;0.988300000000000;0.988090000000000;0.987890000000000;0.987690000000000;0.987500000000000;0.987310000000000;0.987120000000000;0.986940000000000;0.986760000000000;0.986580000000000;0.986410000000000;0.986240000000000;0.986080000000000;0.985920000000000;0.985770000000000;0.985620000000000;0.985470000000000;0.985330000000000;0.985190000000000;0.985060000000000;0.984930000000000;0.984810000000000;0.984690000000000;0.984570000000000;0.984460000000000;0.984360000000000;0.984260000000000;0.984160000000000;0.984070000000000;0.983990000000000;0.983910000000000;0.983830000000000;0.983760000000000;0.983700000000000;0.983630000000000;0.983580000000000;0.983530000000000;0.983480000000000;0.983440000000000;0.983400000000000;0.983370000000000;0.983350000000000;0.983330000000000;0.983310000000000;];

% allocating the variables
meanReflectance = zeros(numberOfImages,242);
stdReflectance = zeros(numberOfImages,242);
sunElevation = zeros(numberOfImages,1);
sunZenith = zeros(numberOfImages,1);
sunAzimuth = zeros(numberOfImages,1);
sensorLookAngle = zeros(numberOfImages,1);
acquisitionDate = cell(numberOfImages,1);
dateOfYear = zeros(numberOfImages,1);
decimalYear = zeros(numberOfImages,2);
% store empty acquisition for removal
empty_acquisition = [];
%%
for dateAcqui = 1:numberOfImages
     imageLocation = fullfile(dir_list(dateAcqui).folder,dir_list(dateAcqui).name,product);
     % print to terminal to see the processing
     warning('Directory Processing : %s',imageLocation);
     
     if exist(imageLocation,'dir')~= 7
       warning('on')
       warningText = strcat('Folder not Found', imageLocation);
       warning(warningText);
       %emptyDir(dateAqui) = 1;
       empty_acquisition = [empty_acquisition, dateAcqui];
       continue;
     end
     
    %searching the directory for the files
    file_list = ls(imageLocation);
    name = file_list(3,:);
    
    %reading the common string in the filename
    baseName = name(1:22);
    
    %clearing the file_list
    clear file_list
    
   %reading the mtl file using the built in function name MTL_parser_L8
    mtlFileName = fullfile(imageLocation,strcat(baseName,'_','MTL_',product,'.txt'));
    if exist(mtlFileName,'file') ~=2
        warningText = strcat('MTL file not Found ',mtlFileName);
        warning(warningText);
        empty_acquisition = [empty_acquisition, dateAcqui];
        continue;
    end
    MTLList = MTLParserHyperion(mtlFileName);
    
    % product information extract
    sunElevation(dateAcqui) = double(MTLList.PRODUCT_PARAMETERS.SUN_ELEVATION);
    sunZenith(dateAcqui) = (90 - sunElevation(dateAcqui));
    sunAzimuth(dateAcqui) = double(MTLList.PRODUCT_PARAMETERS.SUN_AZIMUTH);
    sensorLookAngle(dateAcqui) = double( MTLList.PRODUCT_PARAMETERS.SENSOR_LOOK_ANGLE);
    acquisitionDate{dateAcqui} = MTLList.PRODUCT_METADATA.ACQUISITION_DATE;
    [doy, dp] = date2doy(datenum(acquisitionDate{dateAcqui}));
    dateOfYear(dateAcqui) = doy;
    decimalYear(dateAcqui,1) = dp;
    decimalYear(dateAcqui,2) = str2double(erase(MTLList.PRODUCT_METADATA.ACQUISITION_DATE,'-'));
% %     ULMAPX = MTLList.PRODUCT_METADATA.PRODUCT_UL_CORNER_MAPX;
% %     URMAPX = MTLList.PRODUCT_METADATA.PRODUCT_UR_CORNER_MAPX;
% %     ULMAPY = MTLList.PRODUCT_METADATA.PRODUCT_UL_CORNER_MAPY;
% %     LRMAPY = MTLList.PRODUCT_METADATA.PRODUCT_LR_CORNER_MAPY;
% %     lines = MTLList.PRODUCT_METADATA.PRODUCT_LINES;
% %     samples = MTLList.PRODUCT_METADATA.PRODUCT_SAMPLES;
% %     startTime = MTLList.PRODUCT_METADATA.START_TIME;
% %     endTime = MTLList.PRODUCT_METADATA.END_TIME;
    
    %scaling factors for VNIR and SWIR regions
    scalingFactorVNIR = 40;
    scalingFactorSWIR = 80;
    
    % Setup ROI region
    % Mapx=[XMAP1;XMAP2;XMAP2;XMAP1;XMAP1];
    % Mapy=[YMAP1;YMAP1;YMAP2;YMAP2;YMAP1];
    x_vec = [UL_x; UR_x; LR_x; LL_x; UL_x];
    y_vec = [UL_y; UR_y; LR_y; LL_y; UL_y];
    
    
    %VNIR Bands
    for  k = 8:55
        % read the Hyperion Images
        imageName  = fullfile(imageLocation,strcat(baseName,'_','B',num2str(k,'%03d'),'_',product,'.TIF'));
        if(exist(imageName, 'file') == 0)
            imageName  = fullfile(imageLocation,strcat(baseName,'_','B',num2str(k,'%03d'),'_','_L1GST.TIF'));
        end
        % reading the image 
        [image,R] = geotiffread(imageName);
        
        % CONVERTING THE MAP CO-ORDINATES TO THE pixel
       % [imageRow, imageCol] = map2pix(R, Mapx, Mapy);
        [imageRow, imageCol] = map2pix(R, x_vec, y_vec);
        % PROVIDES THE GEO INFORMATION ABOUT iMAGES
        imageInfo = geotiffinfo(imageName);
        
        % CREATING THE MASK OF THE rEGION OF iNTEREST
        region = poly2mask(imageCol,imageRow,imageInfo.Height,imageInfo.Width);
        
        %APPLYING THE MASK TO THE iMAGE
        maskedImage = double(image).*double(region);
        
        %TAKING THE MEAN OF THE MASKED REGION
        imageMean = nanmean(nonzeros(maskedImage(:)));
        
        %TAKING THE STANDARD DEVIATION OF THE MASKED REGION
        imageStd = nanstd(nonzeros(maskedImage(:)));
        
        % MEAN REFLECTANCE OF THE IMGAE
        % https://eo1.usgs.gov/faq/question?id=21
%         meanReflectance(dateAcqui,k) = (imageMean/scalingFactorVNIR)*pi*sundistance(doy)^2 ...
%             /(esun(k)*sind(sunElevation(dateAcqui))*cosd(sensorLookAngle(dateAcqui)));
         meanReflectance(dateAcqui,k) = (imageMean/scalingFactorVNIR)*pi*sundistance(doy)^2 ...
            /(esun(k)*sind(sunElevation(dateAcqui))*cosd(sensorLookAngle(dateAcqui)));
        
        % STANDARD DEVIATION OF THE IMGAE
%         stdReflectance(dateAcqui,k)= (imageStd/scalingFactorVNIR)*pi*sundistance(doy)^2 ...
%             /(esun(k)*sind(sunElevation(dateAcqui))*cosd(sensorLookAngle(dateAcqui)));
         stdReflectance(dateAcqui,k)= (imageStd/scalingFactorVNIR)*pi*sundistance(doy)^2 ...
            /(esun(k)*sind(sunElevation(dateAcqui))*cosd(sensorLookAngle(dateAcqui)));
    end
    
    %SWIR Bands
    for k = 77:224
         % read the Hyperion Images
        imageName  = fullfile(imageLocation,strcat(baseName,'_','B',num2str(k,'%03d'),'_',product,'.TIF'));
        if(exist(imageName, 'file') == 0)
            imageName  = fullfile(imageLocation,strcat(baseName,'_','B',num2str(k,'%03d'),'_','_L1GST.TIF'));
        end
        [image,R] = geotiffread(imageName);
        [imageRow, imageCol] = map2pix(R, x_vec, y_vec);
        imageInfo = geotiffinfo(imageName);
        %region = poly2mask(imageCol,imageRow,imageInfo.Height,imageInfo.Width);
        maskedImage = double(image).*double(region);
        imageMean = nanmean(nonzeros(maskedImage(:)));
        imageStd = nanstd(nonzeros(maskedImage(:)));
        meanReflectance(dateAcqui,k) = (imageMean/scalingFactorSWIR)*pi*sundistance(doy)^2 ...
            /(esun(k)*sind(sunElevation(dateAcqui))*cosd(sensorLookAngle(dateAcqui)));
        stdReflectance(dateAcqui,k) = (imageStd/scalingFactorSWIR)*pi*sundistance(doy)^2 ...
            /(esun(k)*sind(sunElevation(dateAcqui))*cosd(sensorLookAngle(dateAcqui)));
    end
end
meanReflectance = [transpose(hyperionbandwavelengths);meanReflectance];
stdReflectance  = [transpose(hyperionbandwavelengths);stdReflectance];
% Removing the bands that  are not required 
bandNotrequired  = [1:7,56:76,225:242];
meanReflectance(:,bandNotrequired) = [];
stdReflectance(:,bandNotrequired) = [];

endTime = datetime('now');

timeTaken = endTime-startTime;
% save ('hyperion_over_Libya_4_view_solar.mat','meanReflectance','stdReflectance', ...
%     'empty_acquisition','sunElevation','sunZenith','sunAzimuth',...
%     'sensorLookAngle','acquisitionDate','dateOfYear','decimalYear','ROI'); %

clearvars -except meanReflectance
