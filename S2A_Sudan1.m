%%%%% Sudan 1 temporal trend and other analysis;
clc
clear all;
if ispc
    location='Z:\ImageDrive\Sentinel\MSI-A\P177\R045';
elseif isunix
    localtion='~/zdrive/ImageDrive/Sentinel/MSI/P177/R045';
end

% Get a list of all files and folders in this folder.
files = dir(location);
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);

subFolders(1:2) =[];
%subFolders(9)=[];

length=size(subFolders);

%sudan 1
% 
% lat1=[561570 584250 584250 561570 561570];
% lon1=[2400000 2400000 2367450 2367450 2400000];

% 
lat=[499980 609780 609780 499980 499980];
lon=[2400000 2400000 2290200 2290200 2400000];

% % % Harika
lat1=[561570 584250 584250 561570 561570];
lon1=[24005850 24005850 2367450 2367450 24005850];
% % 
% SBAF=[0.995 0.9722 1.010167 0.97837 0.999832 0.99876 1.000109];

Resolution=[60 10 10 10 20 20 20 10 60 60 20 20 20];

Q_value=10000;

Refl_S2A=zeros(length(1),7);

for i=1:1:length(1)
    file_s=dir(fullfile(location,subFolders(i).name));
    file_s(1:2)=[];
    
    Acquisition_dates=subFolders(i).name;
    [doy,fraction] = date2doy(datenum(Acquisition_dates,'yyyymmdd'));
    dsl= datenum(Acquisition_dates, 'yyyymmdd') - datenum(num2str('2015-06-23'), 'yyyy-mm-dd');
    dsl_vect(i)=dsl;
    doy_vect(i)=doy;
    
    j1=0;
    for j=[1 2 3 4 13 11 12]
        j1=j1+1;
        image_file_name= file_s(j).name;
        image_file_location=fullfile(location,subFolders(i).name,image_file_name);
        
        Image=imread(image_file_location);
        
        x1=(lon1(1)-lon(1))/(-Resolution(j));
        y1=(lat1(1)-lat(1))/(Resolution(j));
        x2=(lon1(2)-lon(1))/(-Resolution(j));
        y2=(lat1(2)-lat(1))/(Resolution(j));
        x3=(lon1(3)-lon(1))/(-Resolution(j));
        y3=(lat1(3)-lat(1))/(Resolution(j));
        x4=(lon1(4)-lon(1))/(-Resolution(j));
        y4=(lat1(4)-lat(1))/(Resolution(j));
        X=[x1 x2 x3 x4];
        Y=[y1 y2 y3 y4];
        
        Image_mask=roipoly(Image,Y,X);
        
        if j==2
            my_mask=Image_mask;
        end
        Image=double(Image).*double(Image_mask);
        
        Field_data=sum(Image(:)==0);
        [row,column]=size(Image);
        
        Reflectance=sum(sum(double(Image)))/(row*column-Field_data)/Q_value;
        Refl_S2A(i,j1)=Reflectance;
        STD_Refl_S2A(i,j1)=nanstd(nanstd(Image))/(row*column-Field_data);
        
    end
    
    
    %%%% View and Azimuth angle information collection for the selected ROI;
    Sun_Zenith=double(imread(fullfile(location,subFolders(i).name,'Sun_Zenith.png'))).*double(my_mask);
    SZA_ROI_S2A(i)=mean2(Sun_Zenith(Sun_Zenith~=0));
    
    Sun_Azimuth=double(imread(fullfile(location,subFolders(i).name,'Sun_Azimuth.png'))).*double(my_mask);
    SAA_ROI_S2A(i)=mean2(Sun_Azimuth(Sun_Azimuth~=0));
    
    View_Zenith=double(imread(fullfile(location,subFolders(i).name,'View_Zenith.png'))).*double(my_mask);
    VZA_ROI_S2A(i)=mean2(View_Zenith(View_Zenith~=0));
    
    View_Azimuth=double(imread(fullfile(location,subFolders(i).name,'View_Azimuth.png'))).*double(my_mask);
    VAA_ROI_S2A(i)=mean2(View_Azimuth(View_Azimuth~=0));
end

%Decimal year estimation from acquisition dates;
str={subFolders.name};
str=str';
sdn = datenum( str, 'yyyymmdd' );
date=datestr( sdn, 'yyyy-mm-dd' );
DateString = date;
formatIn = 'yyyy-mm-dd';
DateVector = datevec(DateString,formatIn);
doy_vect=date2doy(sdn);
Year_S2A=DateVector(:,1);
Decimal_Year_S2A=Year_S2A+doy_vect/365;


% Find scene number for Cloudy scene exclusion again
x = Decimal_Year_S2A;
dt = datetime(floor(x), 1, 1) + years(x-floor(x));
[Y, M, D] = ymd(dt);

% No outlier has been removed from here as the suspected scenes from 2
% sigma consideration are not cloudy scenes found from visual inspection

reflectance=Refl_S2A;

% temporal uncertaintry without BRDF correction;
for p = 1:7
M(p)= mean(reflectance(:,p));
    U(p)= std(reflectance(:,p))/M(p)*100;
end


%%%%%% BRDF correction of the Mean TOA reflectance using 4 angle linear BRDF regression model;
y1=sind(SZA_ROI_S2A).*sind(SAA_ROI_S2A);
x1=sind(SZA_ROI_S2A).*cosd(SAA_ROI_S2A);
y2=sind(VZA_ROI_S2A).*sind(VAA_ROI_S2A);
x2=sind(VZA_ROI_S2A).*cosd(VAA_ROI_S2A);

n=size(reflectance);

for i=1:1:7
    tbl=table(reflectance(:,i),y1,x1,y2,x2,'VAriableNAmes',{'R' 'u' 'v' 'w' 'x'});
    
    lm=fitlm(tbl,'R~u+v+w+x');
    
    % sse(i)=lm.SSE;
    
    b=lm.Coefficients.Estimate;
    
    % reference angles chosen from the average of angles of S2A
    y11=sind(35).*sind(120);
    x11=sind(35).*cosd(120);
    y12=sind(4).*sind(100);
    x12=sind(4).*cosd(100);
    
    
    ref_reflc(i)= b(1)+b(2)*y11+b(3)*x11+b(4)*y12+b(5)*x12; % reference reflectance value;
    
    for j = 1:1:n
        predicted_reflc(j,i)= b(1)+b(2)*y1(j)+b(3)*x1(j)+b(4)*y2(j)+b(5)*x2(j); % predicted reflectnace value:
        
        corrected_reflc(j,i) = reflectance(j,i)*ref_reflc(i)/predicted_reflc(j,i); % corrected reflectance value;
        
    end
end

% BRDF uncertainty calculation (percentage)
S2A_brdfUncertainty = abs(((Refl_S2A-predicted_reflc)./Refl_S2A)*100);

% temporal uncertaintry without BRDF correction;
for p = 1:7
Mc(p)= mean(corrected_reflc(:,p));
    Uc(p)= std(corrected_reflc(:,p))/Mc(p)*100;
end

%%%%%%%%% Plotting temporal trend of TOA reflectance;
COL = {'b','c','g','r','m','[0.6667 0 0]','k'};
col = {'[0.5 0.8 0.2]', '[0.2 0.6 0.8]', '[0.9 0.2 0.6]','[0.3 0.6 0.1]','[0.9 0.5 0.7]','[0.8 0.8 0.2]','[0.4 0.5 0.6]'};
Bands = {'CA','Blue','Green','Red','NIR','SWIR1','SWIR2'};

for i = 1:7
    figure(i)
    
    plot(Decimal_Year_S2A,reflectance(:,i),'LineStyle','none','Color',COL{i},'MarkerSize',20,'Marker','h','MarkerFaceColor',COL{i})
    xlabel('Decimal Year','FontSize',16,'FontWeight','bold','Color','k')
    ylabel('TOA Reflectance','FontSize',16,'FontWeight','bold','Color','k')
    title(['TOA Reflectance vs Decimal Year for ' Bands{i} '_'  'Band'],'FontSize',20);
    ylim ([mean(reflectance(:,i))-0.25   mean(reflectance(:,i))+0.25])
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 30);
    grid on;
    hold on
    plot(Decimal_Year_S2A,corrected_reflc(:,i),'LineStyle','none','Color',col{i},'MarkerSize',20,'Marker','o')
    legend({'Before Correction' 'After Correction' })
    hold off

end   

bandNames = {'CA','Blue','Green','Red','NIR','SWIR1','SWIR2'};

% Shapiro-Wilk test of normality of the lifetime meanTOA reflectance after BRDF using swtest function from file exchange
for band = 1:7
    [H(:,band), pValue(:,band), W(:,band)] =  swtest(corrected_reflc(:,band),0.05); % Shapiro wilk test
    h(:,band) = adtest(corrected_reflc(:,band)); %%%%%% Anderson-Darling test
    %%% Linear relationship check between decimal year and TOA reflectance
    correlation_matrix = corrcoef(Decimal_Year_S2A(:,1),corrected_reflc(:,band));
    % plot of normality check by histogram
    figure (band)
    histogram (corrected_reflc(:,band),20,'Edgecolor','black','LineWidth',1)
    xlabel('Value of TOA Reflectance', 'FontSize',20);
    ylabel('Frequency of Occurrence for TOA Reflectance','FontSize',20);
    title(['Histogram of TOA Reflectance' '-' bandNames{band} '-' 'Band'],'FontSize',20);
    set(gca, 'FontSize', 20)
end
