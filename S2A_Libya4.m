%%%%% Libya 4 temporal trend and other analysis;
clc
clear all;
if ispc
    location='Z:\ImageDrive\Sentinel\MSI\P181\R040';
elseif isunix
    localtion='~/zdrive/ImageDrive/Sentinel/MSI/P181/R040';
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

%Libya 4;
lat1=[723825 743355 743805 724245 723825];
lon1=[3171375 3171825 3149685 3149325 3171375];
% % % %
lat=[699960 809760 809760 699960 699960];
lon=[3200040 3200040 3090240 3090240 3200040];


Resolution=[60 10 10 10 20 20 20 10 60 60 20 20 20];

Q_value=10000;

Refl=zeros(length(1),7);

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
        
        Field_data=sum(Image(:)==0); %%he meant filled data...holy crap. guys please at least use good names
        [row,column]=size(Image);
        
        Reflactance=sum(sum(double(Image)))/(row*column-Field_data)/Q_value;
        Refl(i,j1)=Reflactance;
        STD_Refl(i,j1)=nanstd(nanstd(Image))/(row*column-Field_data);
        
    end
    %%%% View and Azimuth angle information collection for the selected ROI;
    Sun_Zenith=double(imread(fullfile(location,subFolders(i).name,'Sun_Zenith.png'))).*double(my_mask);
    SZA_ROI(i)=mean2(Sun_Zenith(Sun_Zenith~=0));
    
    Sun_Azimuth=double(imread(fullfile(location,subFolders(i).name,'Sun_Azimuth.png'))).*double(my_mask);
    SAA_ROI(i)=mean2(Sun_Azimuth(Sun_Azimuth~=0));
    
    View_Zenith=double(imread(fullfile(location,subFolders(i).name,'View_Zenith.png'))).*double(my_mask);
    VZA_ROI(i)=mean2(View_Zenith(View_Zenith~=0));
    
    View_Azimuth=double(imread(fullfile(location,subFolders(i).name,'View_Azimuth.png'))).*double(my_mask);
    VAA_ROI(i)=mean2(View_Azimuth(View_Azimuth~=0));
end

