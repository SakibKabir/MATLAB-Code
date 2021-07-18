clc
clear all;



%cd('Y:\Libya 1')
if ispc
    location='Z:\ImageDrive\Sentinel\MSI\P043\R033';
elseif isunix
    localtion='~/zdrive/ImageDrive/Sentinel/MSI/P187/R043';
end

% Get a list of all files and folders in this folder.
files = dir(location);
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags);

subFolders(1:2) =[];
%subFolders(9)=[];

length=size(subFolders)



% %%%Site Selection and ROI
% SS=input('Input the first letter of site and then corresponding number \n for Libya 1, input L1 \n Site=')
% Site={'L4' 'L1' 'E1' 'N1' 'N2' 'S1'}
% Sites={'Libya4' 'Libya1' 'Egypt1' 'Niger1' 'Niger2' 'Sudan1'}
% x = strmatch(SS, Site)
% 
% fileID = fopen(['Z:\PNP\ROI_Location\' Sites{x} '_ROI_OPT.txt'],'r');
% % 
% 
% C = textscan(fileID,'%d %d' ,'CommentStyle','//');
% fclose(fileID);
% 
% C=cell2mat(C)




%Extent for ROI SDSU 
% lat=[600000 709800 709800 600000 600000];
% lon=[ 5000040  5000040 4890240 4890240  5000040];
% 
% lat1=[678200 678440 678440 678200 678200];
% lon1=[4906740 4906740 4906620 4906620 4906740];




%%Crater Lake
% lat1=[571695 574965 574965 571695 571695];
% lon1=[4756635 4756635 4753245 4753245 4756635];
% %
% lat=[499980 609880 609880 499980 499980];
% lon=[4800000  4800000 4690200 4690200  4800000];
% 



%%Lake Tahoe
lat1=[751125 762315  762315  751125 751125];
lon1=[4337625 4337625 4326405 4326405 4337625];
% 
lat=[699960 809760 809760 699960 699960];
lon=[4400040 4400040 4290240 4290240 4400040];

SBAF=[1.02 1.08 0.982 1.018 1.005 0.998 0.998];


%%Libya 1
% lat1=[330150 365070 365070 330150 330150];
% lon1=[2750850 2750850 2716860 2716860 2750850];
% 
% lat=[300000 409800 409800 300000 300000];
% lon=[2800020 2800020 2690220 2690220 2800020];
% 
% SBAF=[1.013 0.967 1.018 0.977 0.999 0.998 0.999];


% % %Libya 4
%1
% lat1=[731565 797445 797445 731565 731565];
% lon1=[3189885 3189885 3143865 3143865 3189885];
% 2 (Harika)
% lat1=[723825 743355 743805 724245];
% lon1=[3171375 3171825 3149685 3149325];
% % % % % 
% lat=[699960 809760 809760 699960 699960];
% lon=[3200040 3200040 3090240 3090240 3200040];
% 
% SBAF=[1.008 0.9723 1.005 0.978 0.999 0.998 0.998]; %%changed blue band
% sbaf

% % % % %%Niger 2
% lat1=[644190 677670 677670 644190 644190];
% lon1=[2375910 2375910 2350590 2350590 2375910];
% 
% lat=[600000 709800 709800 600000 600000];
% lon=[2400000 2400000 2290200 2290200 2400000];
% SBAF=[0.994 0.9727 1.010 0.9791 1.0002 0.99915 1.00005113];


% %sudan 1
% % 
% lat1=[561570 584250 584250 561570 561570];
% lon1=[2400000 2400000 2367450 2367450 2400000];
% 
% % 
% lat=[499980 609780 609780 499980 499980];
% lon=[2400000 2400000 2290200 2290200 2400000];
% % % %%Harika
% lat1=[561570 584250 584250 561570 561570];
% lon1=[24005850 24005850 2367450 2367450 24005850];
% % 
% SBAF=[0.995 0.9722 1.010167 0.97837 0.999832 0.99876 1.000109];



% %Niger 1
% lat1=[520470 555120 555120 520470  520470 ];
% lon1=[2271210 2271210 2242890 2242890 2271210];
% 
% lat=[499980 609780 609780 499980 499980];
% lon=[2300040 2300040 2190240 2190240 2300040];
% SBAF=[1.013 0.9591 1.01258 0.9787 0.98213 1.0212 0.99828];


% % %Egypt 1
% % % 
% lat1=[431790 462960 462960 431790 431790];
% lon1=[3000930 3000930 2977110 2977110 3000930];
% 
% 
% lat=[399960 509760 509760 399960 399960];
% lon=[3000000 3000000 2890200 2890200 3000000];
% SBAF=[1.012 0.969 1.006 0.9781 0.9994 0.9987 0.9987];




% % %Extent for ROI Volcanic Site
% lat=[699960 809760 809760 699960 699960];
% lon=[ 2800020  2800020 2690220 2690220  2800020];
% 
% lat1=[775652 784770 784770 775652 775652];
% lon1=[2755834 2755834 2749872 2749872 2755834];



Resolution=[60 10 10 10 20 20 20 10 60 60 20 20 20];

Q_value=10000;

 Refl=zeros(length(1),7);
 
 
 
 
 
 for i=1:1:length(1)
     
     
    file_s=dir(fullfile(location,subFolders(i).name));
    file_s(1:2)=[];
     
     %%%%Metadata for DSL
     
 META=xml2struct(fullfile(location,subFolders(i).name,'metadata.xml')); 
 
 pp=size(META.Children(4).Children(4).Children);
 
 SZA(i)=str2num(META.Children(4).Children(4).Children(4).Children(2).Children.Data);
                 
 SAA(i)=str2num(META.Children(4).Children(4).Children(4).Children(4).Children.Data );
             
 
 
  time_datenum=datenum(subFolders(i).name,'yyyymmdd');
     dol_datenum=datenum('2015-06-23','yyyy-mm-dd');
     dol_datenum1=datenum('2013-02-11','yyyy-mm-dd');
     
     DSLS2(i)=time_datenum-dol_datenum;
 DSLwr2l8(i)=time_datenum-dol_datenum1;
 
 
 
 i1=0;
% % % % % % % % for ii=[2 8 10 12 26 6 24]
% % % % % % % %     
% % % % % % % %     i1=i1+1;
% % % % % % % %     VZA(i,i1)=str2num(META.Children(4).Children(4).Children(pp(2)-1).Children(ii).Children(2).Children.Data);
% % % % % % % %     VAA(i,i1)=str2num(META.Children(4).Children(4).Children(pp(2)-1).Children(ii).Children(4).Children.Data);
% % % % % % % %     
% % % % % % % %     
% % % % % % % % end
 
 
 
 
 
 
 j1=0;
 for j=[1 2 3 4 13 11 12]
       j1=j1+1; 
      image_file_name= file_s(j).name;
        image_file_location=fullfile(location,subFolders(i).name,image_file_name)
        
        
        
        
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
    
%  
% xmin=(lon1(1)-lon(1))/(-Resolution(j));
% ymin=(lat1(1)-lat(1))/(Resolution(j));
% width=(lat1(2)-lat1(1))/Resolution(j)-1;
% height=(lon1(3)-lon1(1))/(-Resolution(j))-1;
%         
% Image=imcrop(Image,[xmin ymin width height]);
        

        Field_data=sum(Image(:)==0);
[row,column]=size(Image);


Reflactance=sum(sum(double(Image)))/(row*column-Field_data)/Q_value;


        Refl(i,j1)=Reflactance;
        
        
 end
     
 
 Sun_Zenith=double(imread(fullfile(location,subFolders(i).name,'Sun_Zenith.png'))).*double(my_mask);
 SZA_ROI(i)=mean2(Sun_Zenith(Sun_Zenith~=0));
 
 Sun_Azimuth=double(imread(fullfile(location,subFolders(i).name,'Sun_Azimuth.png'))).*double(my_mask);
 SAA_ROI(i)=mean2(Sun_Azimuth(Sun_Azimuth~=0));
 
 View_Zenith=double(imread(fullfile(location,subFolders(i).name,'View_Zenith.png'))).*double(my_mask);
 VZA_ROI(i)=mean2(View_Zenith(View_Zenith~=0));
 
 View_Azimuth=double(imread(fullfile(location,subFolders(i).name,'View_Azimuth.png'))).*double(my_mask);
 VAA_ROI(i)=mean2(View_Azimuth(View_Azimuth~=0));
 

 
 end
 
 

 
  PPPPP=Refl
 %%%%SBAF correction
 for i=1:1:7
     
     
 Refl(:,i)= Refl(:,i)*SBAF(i);
     
     
     
 end
 
 
 
 
y1=sind(SZA_ROI').*sind(SAA_ROI');
x1=sind(SZA_ROI').*cosd(SAA_ROI');
 
y2=sind(VZA_ROI').*sind(VAA_ROI');
x2=sind(VZA_ROI').*cosd(VAA_ROI'); 
 
% y2=y2;
% x2=x2;

% Refl=Refl1
 
 Bands={'CA','Blue','Green','Red','NIR','SWIR 1','SWIR 2'};
 
 for i=1:1:7
%%%%
T=table(Refl(:,i),y1,x1,y2,x2,'VariableNames',{'Reflectance' 'y1' 'x1' 'y2' 'x2'}) 
writetable(T,[Bands{i} '.csv'])
%type ['Band' num2str(i) '.csv']

 end
 
 
 
T=table(DSLS2','VariableNames',{'DSL'}) 
writetable(T,'DSL.csv')








 
%  %%%%Deletion
% % % % % l=0
 for i=[14]

Refl(i,:)=[];
% %  Refl1(i,:)=[];
%  SZA(i)=[];
%  SAA(i)=[];
%  VZA(i,:)=[];
%  VAA(i,:)=[];
%  DSL(i)=[];

 SZA_ROI(i)=[];
 SAA_ROI(i)=[];
 VZA_ROI(i)=[];
 VAA_ROI(i)=[];
 DSLS2(i)=[];
 DSLwr2l8(i)=[];
% l=l+1;
 end
% % % % % % % %  
%  

 
 
 
%%%%Saving data
T=table(DSLS2',Refl(:,1),Refl(:,2),Refl(:,3),Refl(:,4),Refl(:,5),Refl(:,6),Refl(:,7),'VariableNames',{'DSL' 'CA' 'Blue' 'Green' 'Red' 'NIR' 'SWIR_1' 'SWIR_2'})
 
 %writetable(T,'L4.xls')
 writetable(T,'L4_S2.csv')
 
 clear all
 
 
 
 
 
 
 
 
 SBAF=[0.9995 0.9723 1.0031 0.9753 1.0002 0.9985 0.9984;0.9994 0.9670 1.0197 0.9776 0.9993 0.9987 0.9997;...
    0.9995 0.9690 1.0046 0.9792 0.9984 0.9989 0.9981;0.9995 0.9591 1.0113 0.9780 0.9991 0.9988 0.9995;...
    0.995 0.9727 1.0104 0.9804 1.0004 0.9956 0.9967;0.9995 0.9722 1.0082 0.9787 0.999 0.9989 1.0018];
 
 SBAF=SBAF(4,:)

 
 PPPPP=Refl
 %%%%SBAF correction
 for i=1:1:7
     
     
 Refl(:,i)= Refl(:,i)*SBAF(i);
     
     
     
 end
 

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 %positive and negative sza for modeling
 
 [~,b]=size(SZA)
 
 %Ref_VAA=mean2(VAA(:,1));
 
 
 for j=1:1:7
     
      Ref_VAA=mean2(VAA(:,j));
      
 for i=1:1:b
     
 if VAA(i,j)>=Ref_VAA  
     
   VZA(i,j)=-VZA(i,j); 
 end
     
 end
 
 end
 
 
 
 
 
 
 
 
 
 
 
 
 
 %%%% Polar plot of the sun and satellite geometry
 polarplot(deg2rad(SAA),SZA,'bp', 'MarkerSize',12)
 hold on 
 polarplot(deg2rad(VAA(:,1)),VZA(:,1),'bo')
 hold on
 pax = gca;
pax.ThetaAxisUnits = 'degree';
pax.FontSize = 18;
 title('Sun and Sensor geometry for PICS')
 clear all
 
%  xlabel('Sun Zenith Angle','FontSize',16,'FontWeight','bold','Color','k')
%  ylabel('TOA Reflectance','FontSize',16,'FontWeight','bold','Color','k')
 title('Sun Zenith angle and View Zenith angles of different sites','FontSize',18);
 
 
 
 COL={'b','c','g','r','m','[0.6667 0 0]','k'};
 
 for i=1:7
     %c=[COL{i} '*']
 polarplot(deg2rad(SZA),Refl1(:,i),'LineStyle','none','Marker','o','markerfacecolor',COL{i},'markeredgecolor',COL{i})
 hold on 
 pax = gca;
pax.ThetaAxisUnits = 'degree';
 %set(L,'markerfacecolor',COL{i})
 pax.FontSize = 12;
 end
 legend(Bands,'Location','northeast')
 title('Reflectance vs SZA for Libya 4')
 
 clear all
 
 
 
 
 
 
 
 %%%%%Writing in xl file without ROI
 
y1=sind(SZA').*sind(SAA');
x1=sind(SZA').*cosd(SAA');
 
y2=sind(VZA').*sind(VAA');
x2=sind(VZA').*cosd(VAA'); 
 
y2=y2';
x2=x2';

%  Refl=Refl1
 
 Bands={'CA','Blue','Green','Red','NIR','SWIR 1','SWIR 2'};
 
 for i=1:1:7
 %%%%
T=table(Refl1(:,i),y1,x1,y2(:,i),x2(:,i),'VariableNames',{'Reflectance' 'y1' 'x1' 'y2' 'x2'}) 
writetable(T,[Bands{i} '.csv'])
%type ['Band' num2str(i) '.csv']

 end
 
 
 
 
 
  %%%%%Writing in xl file without ROI
  

%  Refl=Refl1
 
 Bands={'CA','Blue','Green','Red','NIR','SWIR 1','SWIR 2'};
 
 for i=1:1:7
 %%%%
T=table(Refl(:,i),y1,x1,y2(:,i),x2(:,i),'VariableNames',{'Reflectance' 'y1' 'x1' 'y2' 'x2'}) 
writetable(T,[Bands{i} '.csv'])
%type ['Band' num2str(i) '.csv']

 end
 
 
 
 
 
 %%%with ROI
 
y1=sind(SZA_ROI').*sind(SAA_ROI');
x1=sind(SZA_ROI').*cosd(SAA_ROI');
 
y2=sind(VZA_ROI').*sind(VAA_ROI');
x2=sind(VZA_ROI').*cosd(VAA_ROI'); 
 
% y2=y2;
% x2=x2;

%  Refl=Refl1
 
 Bands={'CA','Blue','Green','Red','NIR','SWIR 1','SWIR 2'};
 
 for i=1:1:7
 %%%%
T=table(Refl1(:,i),y1,x1,y2,x2,'VariableNames',{'Reflectance' 'y1' 'x1' 'y2' 'x2'}) 
writetable(T,[Bands{i} '.csv'])
%type ['Band' num2str(i) '.csv']

 end
 
 
 
  Bands={'CA','Blue','Green','Red','NIR','SWIR 1','SWIR 2'};
 
 for i=1:1:7
 %%%%
T=table(Refl1(:,i),y1,x1,y2,x2,'VariableNames',{'Reflectance' 'y1' 'x1' 'y2' 'x2'}) 
writetable(T,[Bands{i} '.csv'])
%type ['Band' num2str(i) '.csv']

 end
 
 
 
 
 
 
 % Plotting new corrected 
 
 
 
 R1=Refl;
DSL1=DSLS2;

[L ~]=size(R1);
%clearvars -except R1 DSL1

F=[.4241 2.095e-2 -3.725e-5 -3.885 -4.486;0.470104 0.016789 0.005526 -4.299233 -4.821389;0.505194 -0.002479 0.014410 -3.205847 -3.2589;...
    0.594742 -0.020490 0.027655 -2.016251 -1.660704;0.693943 -0.035774 0.040219 -1.590773 -1.614761;...
    0.634722 -0.084923 0.062618 1.856726 1.110512;0.718941 0.021054 0.051981 -2.174776 -0.370473];

%ref=F(1)+F(2)*(sind(30).*sind(100))+F(3)*(sind(30).*cosd(100))+F(4)*(sind(10).*sind(95))+F(5)*(sind(10).*cosd(95));


%ref=median(R1);
for i=1:1:L%-2

 for j=1:1:7

ref=F(j,1)+F(j,2)*(sind(20).*sind(100))+F(j,3)*(sind(20).*cosd(100))+F(j,4)*(sind(10).*sind(120))+F(j,5)*(sind(10).*cosd(120));
R2(i,j)=(Refl(i,j))*ref/(F(j,1)+F(j,2)*y1(i)+F(j,3)*x1(i)+F(j,4)*y2(i)+F(j,5)*x2(i));



 end
end


mul=[1.0178    0.976    0.9068    0.902    0.98    1.1890    0.7242]

%%Deletion

for i=[42 8]
    
    R1(i,:)=[];
    R2(i,:)=[];
    DSL1(i)=[];
end

Bands={'CA','Blue','Green','Red','NIR','SWIR 1','SWIR 2'};
COL={'b','c','g','r','m','[0.6667 0 0]','k','[0.5 0.75 0.55]'};

for i=1:1:7
figure(i)
plot(DSL1,R1(:,i),'LineStyle','none','Color',COL{i},'Marker','o','MarkerFaceColor',COL{i},'MarkerSize',10); 
 
  mu = mean(R1(:,i));
    st=std(R1(:,i))
% hline = refline([0 mu]);
% hline.Color = 'b';



hold on 
R2(:,i)=R2(:,i)/mul(i);
plot(DSL1,R2(:,i),'LineStyle','none','Color',COL{8},'Marker','^','MarkerFaceColor',COL{8},'MarkerSize',10); 


  mu1 = mean(R2(:,i));
        st1=std(R2(:,i))
% hline = refline([0 mu1]);
% hline.Color = 'r';

% hline = refline([0 ref]);
% hline.Color = 'm';

    xlabel('DSL','FontSize',16,'FontWeight','bold','Color','k')
    ylabel('TOA Reflectance','FontSize',16,'FontWeight','bold','Color','k')
    title(['TOA Reflectance vs DSL for ' Bands{i} ' ' 'Band'],'FontSize',20);
    
  
    
    dim = [.2 .75 .1 .1];
    str = {strcat('\color{black}', '\rho_{Before Correction}=', sprintf('%.3f',mu),'±', sprintf('%.3f', st)),...
        strcat('\color{black}', '\rho_{After Correction}=', sprintf('%.3f',mu1),'±', sprintf('%.3f', st1))};
%annotation('textbox',dim,'String',str)

annotation('textbox',...
    dim,...
    'String',str,...
    'FontSize',20,...
    'EdgeColor','k')

    
    
ylim([mu-0.15 mu+0.15]) 


xlabel('DSL','FontSize',16,'FontWeight','bold','Color','k')
    ylabel('Reflectance','FontSize',16,'FontWeight','bold','Color','k')
    title(['Before and after correction Reflectance of ' Bands{i}  ' Band'],'FontSize',30);
    
    
    
grid on
legend({'Before Correction' 'After Correction'})
   xt = get(gca, 'XTick');
set(gca, 'FontSize', 30)
%aaaaaa(i)=mu1/mu

end





%%Deletion

for i=[42 8]
    
    R1(i,:)=[];
    R2(i,:)=[];
    DSL1(i)=[];
end

 
T=table(DSLS2',Refl(1:7,1),Refl(1:7,2),Refl(1:7,3),Refl(1:7,4),Refl(1:7,5),Refl(1:7,6),Refl(1:7,7),SZA_ROI',SAA_ROI',VZA_ROI',VAA_ROI','VariableNames',{'DSL','CA','Blue','Green','Red','NIR','SWIR_1','SWIR_2', 'SZA', 'SAA', 'VZA', 'VAA'})
writetable(T,['L4_22_s2' '.csv']) 
 
% filename = 'BRDF1.mat';
% save(filename)

 
 
 
 

cd('Y:\Libya 1\BRDF\New\S2')

% Avg_Reflectance=Reflectance;

%Refl=cell2mat(Avg_Reflectance);

SZA=SZA';









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Models%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 
 COL={'b','c','g','r','m','[0.63 0.63 1]','k'};
%  
%  for i=1:1:7
%    %figure(i)
%     plot(SAA,Refl(:,i),'LineStyle','none','Color',COL{i},'Marker','o','MarkerFaceColor',COL{i}); 
%     hold on 
%     xlabel('Sun Azimuth Angle','FontSize',16,'FontWeight','bold','Color','k')
%     ylabel('Reflectance','FontSize',16,'FontWeight','bold','Color','k')
%     title('Reflectance vs Sun Azimuth Angle','FontSize',25);
%      
%  end
% %  saveas(gcf,'SAA') 
% %  saveas(gcf,strcat('SAA','.jpg')); 
%   hold off
%    
 
            Bands={'CA','Blue','Green','Red','NIR','SWIR 1','SWIR 2'};
            COL={'b','c','g','r','m','[0.63 0.63 1]','k'};
%%%%%SZAAA

   for i=1:1:7
   %figure(i)
    plot(DSLS2',Refl(:,i),'LineStyle','none','Color',COL{i},'Marker','o','MarkerFaceColor',COL{i},'MarkerSize',20); 
    hold on 
    xlabel('DSL','FontSize',16,'FontWeight','bold','Color','k')
    ylabel('TOA Reflectance','FontSize',16,'FontWeight','bold','Color','k')
    title(['TOA Reflectance vs DSL for Crater Lake'],'FontSize',18);
    grid on
   end
  xt = get(gca, 'XTick');
set(gca, 'FontSize', 20)
legend(Bands)
 saveas(gcf,'SZA') 
    saveas(gcf,strcat('SZA','.jpg')); 
  hold off
  
  











%%%%%SZAAA

   for i=1:1:7
   %figure(i)
    plot(SZA,Refl1(:,i),'LineStyle','none','Color',COL{i},'Marker','o','MarkerFaceColor',COL{i}); 
    hold on 
    xlabel('Sun Zenith Angle','FontSize',16,'FontWeight','bold','Color','k')
    ylabel('TOA Reflectance','FontSize',16,'FontWeight','bold','Color','k')
    title('TOA Reflectance vs Sun Zenith Angle','FontSize',18);
    grid on
   end
 saveas(gcf,'SZA') 
    saveas(gcf,strcat('SZA','.jpg')); 
  hold off

  
  
  
  % %%%%%VZAA
  
  
  for i=1:1:7
    %figure(i)
    plot((VZA(:,i)),Refl1(:,i),'LineStyle','none','Color',COL{i},'Marker','o','MarkerFaceColor',COL{i}); 
    hold on 
    xlabel('VIew Zenith Angle','FontSize',16,'FontWeight','bold','Color','k')
    ylabel('TOA Reflectance','FontSize',16,'FontWeight','bold','Color','k')
    title('TOA Reflectance vs View Zenith Angle','FontSize',18);
    grid on
  end
    saveas(gcf,'VZA') 
    saveas(gcf,strcat('VZA','.jpg')); 
    hold off
  
  
  %%%%%%VAA
  
  
  for i=1:1:7
    %figure(i)
    plot(VAA(:,i),Refl1(:,i),'LineStyle','none','Color',COL{i},'Marker','o','MarkerFaceColor',COL{i}); 
    hold on 
    xlabel('VIew Azimuth Angle','FontSize',16,'FontWeight','bold','Color','k')
    ylabel('TOA Reflectance','FontSize',16,'FontWeight','bold','Color','k')
    title('TOA Reflectance vs VIew Azimuth Angle','FontSize',18);
    %ylim([limit1(i) limit2(i)]);
    grid on
  end
    saveas(gcf,'VAA') 
    saveas(gcf,strcat('VAA','.jpg')); 
    hold off
  
  
  
  
  %%%%SAA
  for i=1:1:7
    %figure(i)
    plot(SAA,Refl(:,i),'LineStyle','none','Color',COL{i},'Marker','o','MarkerFaceColor',COL{i}); 
    hold on 
    xlabel('Sun Azimuth Angle','FontSize',16,'FontWeight','bold','Color','k')
    ylabel('TOA Reflectance','FontSize',16,'FontWeight','bold','Color','k')
    title('TOA Reflectance vs Sun Azimuth Angle','FontSize',18);
    
  end
   grid on;
   saveas(gcf,'SAA') 
   saveas(gcf,strcat('SAA','.jpg')); 
   hold off


%%%%Deletion
% % % % l=0
 for i=[59 13 6]

 Refl(i,:)=[];
 %Refl1(i,:)=[];
 SZA(i)=[];
 SAA(i)=[];
 VZA(i,:)=[];
 VAA(i,:)=[];
 
 SZA_ROI(i)=[];
 SAA_ROI(i)=[];
 VZA_ROI(i)=[];
 VAA_ROI(i)=[];
 
 
 %DSL(i)=[];
 DSLS2(i)=[];
DSLwr2l8(i)=[];
% l=l+1;
 end
% % % % % % %  
%  

 

 
 
%  
%  
%  
%  for i=1:1:7
% %    figure(i)
% %     plot(SAA,Refl(:,i),'LineStyle','none','Color',COL{i},'Marker','o','MarkerFaceColor',COL{i}); 
% %     %hold on 
% %     xlabel('Sun Azimuth Angle','FontSize',16,'FontWeight','bold','Color','k')
% %     ylabel('Reflectance','FontSize',16,'FontWeight','bold','Color','k')
% %     title('Reflectance vs Sun Azimuth Angle','FontSize',25);
% % %     saveas(gcf,strcat('SAA_B',num2str(i))); 
% %     saveas(gcf,strcat('SAA_B',num2str(i),'.jpg')); 
%     [F, gof] = fit( SAA',Refl(:,i),  'poly2' );
%     
%     %%%Reference reflectance for each band at SAA=90 degree
%     Ref_SAA=90;    %%%reference saa
%     Ref_Refl(i)=F(Ref_SAA);
%     
%     
%     %%%%%%%Correction of reflectance at the reference SZA
%     
%     
%     
%     for j=1:1:length(1)
%         
%         
%         
%         
%         
%       Refl_BRDF_Corrected(j,i)=(Refl(j,i)* Ref_Refl(i))/F(SAA(j)) ;
%         
%         
%         Diff(j,i)=(Refl_BRDF_Corrected(j,i)-Refl(j,i))/Refl(j,i);
%         
%         
%         
%         
%         
%         
%     end
%     
%     
%     %%temporal uncertainty after correction
%       T_uncertainty_SAA(i)=std(Refl_BRDF_Corrected(:,i))/mean(Refl_BRDF_Corrected(:,i));
%     
%     
%     
%     
%     
%     
% %     plot( F,'r-',SAA',Refl(:,i),'bo');
% %     xlabel('Sun Azimuth Angle','FontSize',16,'FontWeight','bold','Color','k')
% %     ylabel('Reflectance','FontSize',16,'FontWeight','bold','Color','k')
% %     title(strcat('Reflectance vs Sun Azimuth Angle of Band',num2str(i)),'FontSize',17);
% %     %%%Linear
% %     saveas(gcf,strcat('SAA_B',num2str(i),'_fit')); 
% %     saveas(gcf,strcat('SAA_B',num2str(i),'_fit','.jpg')); 
%     
%     
%     %quadratic
% %     saveas(gcf,strcat('SAA_B',num2str(i),'_fit_quadratic')); 
% %     saveas(gcf,strcat('SAA_B',num2str(i),'_fit_quadratic','.jpg')); 
% % %     
% %     Slope_SAA(i)=F.p1
% %  R_square_SAA(i)=gof.rsquare
% %   SSE_SAA(i)=gof.sse
%   
%  end
%   hold off
%    
%     
%   
%   
%   
%   
%   
% %%Libya 4 
% limit1=[0.10 0.15 0.20 0.30 0.45 0.55 0.45];
% limit2=[0.35 0.35 0.45 0.55 0.70 0.80 0.75];
% 
% F=[0.1040 1.225e-5 1.174e-3;0.1814 -9.169e-5 5.348e-4;2.499e-1 -2.564e-4 8.099e-4; 0.4103706 -0.0004768 0.0005714;.5754 -7.248e-4 2.646e-4;7.183e-1 -1.14e-3 -3.660e-5;0.4650235 -0.0008991 0.0013627];




%Libya 1 
limit1=[0.17 0.19 0.27 0.40 0.53 0.64 0.54];
limit2=[0.27 0.29 0.37 0.50 0.63 0.74 0.64];

F=[-6.954e-2 7.82e-5 9.847e-4;-1.621e-2 6.296e-5 7.981e-4;8.465e-2 6.56e-5 8.315e-4;4.722e-1 -2.397e-4 1.044e-4;6.144e-1 -5.964e-4 8.452e-5;7.831e-1 -1.096e-3 -1.321e-4;0.5354261 -0.0003894 0.0003111]




%%Egypt 1
% 
% limit1=[0.17 0.19 0.28 0.42 0.52 0.63 0.55];
% limit2=[0.27 0.29 0.38 0.52 0.62 0.73 0.65];
% F=[1.684e-1 9.83e-5 1.876e-4;1.784e-1 -9.237e-5 6.063e-4;2.467e-1 -2.848e-4 6.063e-4;2.467e-1 -2.848e-4 8.727e-4;4.59e-1 -5.292e-4 1.717e-4;0.6514064 -0.0008030 -0.0003747;9.467e-1 -1.289e-3 -2.02e-3;0.7741423 -0.0009804 -0.001331];


%%Niger 1
% limit1=[0.14 0.15 0.26 0.43 0.52 0.63 0.55];
% limit2=[0.24 0.25 0.36 0.53 0.62 0.73 0.65];

% F=[.1142 9.11e-5 4.78e-4;3.879e-1 3.439e-5 -8.879e-4;4.224e-1 1.021e-4 -5.46e-4;.3116 -9.203e-5 8.06e-4;4.204e-1 -4.582e-4 8.541e-4;0.9620505 -0.0012574 -0.0010958;-0.0952049 -0.0006653 0.0036369];
%%Niger 2
% limit1=[0.14 0.15 0.23 0.33 0.45 0.58 0.51];
% limit2=[0.24 0.25 0.33 0.43 0.55 0.68 0.61];
%%%%Golod ase
% F=[1.8041487 0.0000722 -0.0056459;-5.829e-1 2.775e-5 2.723e-3;-2.511 -3.103e-5 9.788e-3;-6.145 -2.842e-4 2.317e-2;4.4665397 -0.0004359 -0.0138798;-4.1369280 -0.0004518 0.0165254];

%%Sudan 1
% limit1=[0.14 0.15 0.23 0.38 0.50 0.63 0.53];
% limit2=[0.24 0.25 0.33 0.48 0.60 0.73 0.63];
% F=[2.448e-1 -7.681e-5 -2.007e-4;.2448 -8.09e-5 -2.078e-4;4.227e-1 -2.307e-5 -6.13e-4;5.954e-1 -1.977e-4 -7.444e-4;.6692 -4.663e-4 -5.405e-4;7.06e-1 -9.008e-4 -3.362e-5;0.5267832 -0.0004807 0.0002696];

%%%%%Algodones dunes
% limit1=[0.16 0.15 0.20 0.30 0.40 0.45 0.4];
% limit2=[0.26 0.25 0.30 0.40 0.50 0.55 0.5];

%%NIger 2
% limit1=[0.18 0.20 0.30 0.40 0.53 0.64 0.65];
% limit2=[0.28 0.30 0.40 0.50 0.63 0.74 0.75];

Bands={'CA','Blue','Green','Red','NIR','SWIR 1','SWIR 2'};


   Refl_BRDF_Corrected=[];
   for i=1:1:7
       
       %%%%Plotting of the normal data without fitting
%    figure(i)
%     plot(SZA,Refl(:,i)','LineStyle','none','Color',COL{i},'Marker','o','MarkerFaceColor',COL{i}); 
%     %hold on 
%     xlabel('Sun Zenith Angle','FontSize',16,'FontWeight','bold','Color','k')
%     ylabel('TOA Reflectance','FontSize',16,'FontWeight','bold','Color','k')
%     title(['TOA Reflectance vs Sun Zenith Angle for',' ',Bands{i}],'FontSize',18);
% ylim([limit1(i) limit2(i)]); 
% %%%saving data
%     saveas(gcf,strcat('SZA_B',num2str(i))) 
% saveas(gcf,strcat('SZA_B',num2str(i),'.jpg')); 




%%fitting of the data
%     [F, gof] = fit( SZA',Refl(:,i),  'poly1' );
    %F=regress(Refl1(:,i),[ones(size(SZA))',SZA',VAA(:,i)],0.05)
    %%Reference reflectance for each band at SZA=30 degree
    Ref_SZA=30;    %%%reference sza
    Ref_Refl(i)=F(i,1)+F(i,2)*10+F(i,3)*250;
    
    Yfit{i}=F(i,1)+F(i,2)*SZA'+F(i,3)*VAA(:,i);
    
    
    %%%%%%Correction of reflectance at the reference SZA
    [X1FIT,X2FIT] = meshgrid(SZA',VAA(:,i));
    Yfit{i}=F(i,1)+F(i,2)*X1FIT+F(i,3)*X2FIT;
    figure(i)
   mesh(X1FIT,X2FIT,cell2mat(Yfit(i))) 
   surf(X1FIT,X2FIT,cell2mat(Yfit(i)))
   hold on;
   plot3(SZA',VAA(:,i),Refl1(:,i),'ro','MarkerFacecolor','r')
    xlabel('Sun Zenith Angle','FontSize',12,'FontWeight','bold','Color','k')
    zlabel('TOA Reflectance','FontSize',12,'FontWeight','bold','Color','k')
ylabel('View Azimuth Angle','FontSize',12,'FontWeight','bold','Color','k')
 title(['BRDF Model for',' ',Bands{i},' ','Band'],'FontSize',18);
%  ylim([limit1(i) limit2(i)]); 
set(get(gca,'xlabel'),'rotation',19);
set(get(gca,'ylabel'),'rotation',-24);
    saveas(gcf,strcat('Surface_model_B',num2str(i))); 
    saveas(gcf,strcat('Surface_model_B',num2str(i),'.jpg')); 
   
   
    [pppp ~]=size(SZA')
    for j=1:1:pppp
        
        
        
        
        
        Refl_BRDF_Corrected(j,i)=(Refl1(j,i)* Ref_Refl(i))/(F(i,1)+F(i,2)*SZA(j)+F(i,3)*VAA(j)) ;
        
        
        Diff(j,i)=(Refl_BRDF_Corrected(j,i)-Refl1(j,i))/Refl1(j,i);
        
        
        
        
        
        
    end
    
    
    %temporal uncertainty after correction
    T_uncertainty(i)=std(Refl_BRDF_Corrected(:,i))/mean(Refl_BRDF_Corrected(:,i));
    
     T_uncertainty_b(i)=std(Refl1(:,i))/mean(Refl1(:,i));
    %%%%%Plotting of the fitted models
    
%     plot( F,'r-',SZA',Refl(:,i),'bo');
%     grid on;
%     xlabel('Sun Zenith Angle','FontSize',16,'FontWeight','bold','Color','k')
%     ylabel('TOA Reflectance','FontSize',16,'FontWeight','bold','Color','k')
%     title(['TOA Reflectance vs Sun Zenith Angle for',' ',Bands{i},' ','Band'],'FontSize',14);
%  ylim([limit1(i) limit2(i)]); 
%     
%     %Linear
%     saveas(gcf,strcat('SZA_B',num2str(i),'_fit')); 
%     saveas(gcf,strcat('SZA_B',num2str(i),'_fit','.jpg')); 
% % %     
%     Slope_SZA(i)=F.p1
%     Bias_SZA(i)=F.p2
%     R_square_SZA(i)=gof.rsquare
%     SSE_SZA(i)=gof.sse

%Quadratic
%  saveas(gcf,strcat('SZA_B',num2str(i),'_fit_quadratic')); 
%     saveas(gcf,strcat('SZA_B',num2str(i),'_fit_quadratic','.jpg')); 
%       Slope_SZA(i)=F.p1
%  R_square_SZA(i)=gof.rsquare
%   SSE_SZA(i)=gof.sse
   end
  hold off
 
  
  

  
  
  
for i=1:1:7

    figure(i)
     plot(DSLS2,Refl(:,i),'bo');
%     
    hold on;
    plot(DSLS2,Refl_BRDF_Corrected(:,i),'ks','MarkerFacecolor','k');
    hold on
%     plot(DSLL8,Ref_cor_l8(:,i),'rd');
    
    grid on;
    
    xlabel('DSL','FontSize',16,'FontWeight','bold','Color','k')
    ylabel('TOA Reflectance','FontSize',16,'FontWeight','bold','Color','k')
    title(['TOA Reflectance vs DSL for',' ',Bands{i}],'FontSize',16);
    ylim([limit1(i)-0.1 limit2(i)+0.1]); 
    
    legend({'S2 Reflectance','BRDF corrected S2 Reflectance'},'Location','northeast')
    
    saveas(gcf,strcat('BRDF_Correcter_B',num2str(i),'_fit')); 
    saveas(gcf,strcat('BRDF_Correcter_B',num2str(i),'_fit','.jpg'));
    
    
end
  
  
  
  
  
  
  
  
  
  
  
  %%Temporal uncertainty
  
%%%%Plot against wavelength

Central_wavelength=[448 482 562 665 865 1610 2200];

% plot(Central_wavelength,Slope_SZA(:),'o')
% ylim([-.003 .0005])

[F2 gof]= fit(Central_wavelength',Slope_SZA(:),'exp2')
  
 plot( F2,'r-',Central_wavelength',Slope_SZA(:),'bo'); 
  grid on
  xlabel('Wavelength in nm','FontSize',16,'FontWeight','bold','Color','k')
    ylabel('Slope coefficient','FontSize',16,'FontWeight','bold','Color','k')
    title(['Slope coefficient vs wavelength'],'FontSize',18);
    ylim([limit1(i) limit2(i)]); 
    
    ylim([-.003 .0005])
    xlim([400 2300])
    
    saveas(gcf,'slope'); 
    saveas(gcf,['slope','.jpg']); 
%     
    plot(Central_wavelength',Slope_SZA(:),'bo'); 
    grid on;
 xlabel('Wavelength in nm','FontSize',16,'FontWeight','bold','Color','k')
    ylabel('Slope coefficient','FontSize',16,'FontWeight','bold','Color','k')
    title(['Slope coefficient vs wavelength'],'FontSize',18);
   % ylim([limit1(i) limit2(i)]); 
    
    ylim([-.003 .0005])
    xlim([400 2300])
    
    saveas(gcf,'slope_without_fit'); 
    saveas(gcf,['slope_without_fit','.jpg']); 
% 
% 
% Refl(23,:)=[];
% SZA(23)=[];

% subFolders(7)=[];
% 
%     

%%%To save the workspace for future use
% % 
% filename = 'BRDF1.mat';
% save(filename)
% 

%%To load teh saved workspace
% 
% load('BRDF.mat')





%%%%cd('Y:\Egypt 1\BRDF\New\S2A')
%%%%%%%%Plot of Temporal uncertainty before and after
 
 c = categorical({'CA','Blue','Green','Red','NIR','SWIR 1','SWIR 2'});
Uncertainty = [T_uncertainty_b(1),T_uncertainty(1);T_uncertainty_b(2),T_uncertainty(2);T_uncertainty_b(3),T_uncertainty(3);T_uncertainty_b(4),T_uncertainty(4);T_uncertainty_b(5),T_uncertainty(5);T_uncertainty_b(6),T_uncertainty(6);T_uncertainty_b(7),T_uncertainty(7)];
bar(c,Uncertainty)
 

   
 xlabel('Bands','FontSize',16,'FontWeight','bold','Color','k')
 ylabel('Temporal Uncertainty','FontSize',16,'FontWeight','bold','Color','k')
 title(['Temporal uncertainty before and after BRDF correction'],'FontSize',15);
 grid on;
 
 
 
 set(gca, 'XTick', [1 2 3 4 5 6 7])
set(gca, 'XTickLabel', {'CA','Blue','Green','Red','NIR','SWIR 1','SWIR 2'},'FontWeight','bold')


 set(gca, 'YTick', [0.005 0.01 0.015 0.02 0.025 0.03 0.035 0.04])
set(gca, 'YTickLabel', {'0.5%','1%','1.5%','2%','2.5%','3%','3.5%','4%'})
legend({'Before Correction','After Correction'},'Location','northwest')

 saveas(gcf,'bar'); 
    saveas(gcf,['bar','.jpg']); 



%  xt = get(gca, 'XTick');
% set(gca, 'FontSize', 12)
%  xlhand = get(gca,'xlabel')
% set(xlhand,'string','X','fontsize',20)
 
%  %%%%Limit
% limit1(3)=.25;
% limit2(3)=.35;
 
 
for i=1:1:6
 
    subplot(2,3,i)
 
    plot(SZA',Refl(:,i),'bo');
    hold on
    plot(SZA',Refl_BRDF_Corrected(:,i),'rd');
    grid on;
    
    
%     
%     v = axis;
% handle=title([Bands{i}])
% set(handle,'Position',[.5 .5]);


  title([Bands{i}],'FontSize',8);
    xlabel('Sun Zenith Angle','FontSize',12,'FontWeight','bold','Color','k')
    ylabel('TOA Reflectance','FontSize',12,'FontWeight','bold','Color','k')
    %legend({'Before','After'},'Location','northwest')
    
    
 ylim([limit1(i) limit2(i)]); 
 
 
 set(gcf,'NextPlot','add');
axes;
h = title('MyTitle');
set(gca,'Visible','off');
set(h,'Visible','on');
 
 
end

% AX=legend('plot 1', 'plot 2','location','northeast','orientation','horizontal');
xlabel('Sun Zenith Angle','FontSize',16,'FontWeight','bold','Color','k')
ylabel('TOA Reflectance','FontSize',16,'FontWeight','bold','Color','k')
title('Before and after BRDF corrected TOA reflectance for different Bands','FontSize',20);






hold off
saveas(gcf,'Comparison'); 
saveas(gcf,['Comparison','.jpg']); 







%%%%%%%for swir 2

i=7;
figure(i+5)
plot(SZA',Refl(:,i),'bo');
hold on
plot(SZA',Refl_BRDF_Corrected(:,i),'rd');
grid on;
    
    
%     
%     v = axis;
% handle=title([Bands{i}])
% set(handle,'Position',[.5 .5]);

title(['Before and after BRDF corrected Reflectance for',' ',Bands{i},' ','Band'],'FontSize',15);
xlabel('Sun Zenith Angle','FontSize',15,'FontWeight','bold','Color','k')
ylabel('TOA Reflectance','FontSize',15,'FontWeight','bold','Color','k')
legend({'Before Correction','After Correction'},'Location','northeast')

    
 ylim([limit1(i) limit2(i)]);  
 
 
 
 saveas(gcf,'Comparison SWIR 2'); 
    saveas(gcf,['Comparison SWIR 2','.jpg']); 
   hold off 
%%Limit
% limit1(7)=.55;
% limit2(7)=.65;
    
    
    
    %%%% Absolute Gain Factor calculation
  for i=1:1:7
      
      
      S2(i)=mean(Refl_BRDF_Corrected(:,i));
      
  end
  
  
  
  
%%%%Saving average reflectance values  
save('S2.mat','S2')

%%%%Saving average reflectance values  
save('Refl_BRDF_Corrected_S2.mat','Refl_BRDF_Corrected')
  
% 
% filename = 'BRDF.mat';
% save(filename)

  
  
  




%%%%%%%%%%%%%%%%%%%%%%%%%%%%2*2 plots



figure(1)
for i=1:1:4
 
    subplot(2,2,i)
 
    plot(SZA',Refl(:,i),'bo');
    hold on
    plot(SZA',Refl_BRDF_Corrected(:,i),'rd');
    grid on;
    
    
%     
%     v = axis;
% handle=title([Bands{i}])
% set(handle,'Position',[.5 .5]);


  title([Bands{i}],'FontSize',8);
    xlabel('Sun Zenith Angle','FontSize',12,'FontWeight','bold','Color','k')
    ylabel('TOA Reflectance','FontSize',12,'FontWeight','bold','Color','k')
    %legend({'Before','After'},'Location','northwest')
    
    
 ylim([limit1(i) limit2(i)]); 
 
 
 set(gcf,'NextPlot','add');
axes;
h = title('MyTitle');
set(gca,'Visible','off');
set(h,'Visible','on');
 
 
end

% AX=legend('plot 1', 'plot 2','location','northeast','orientation','horizontal');
xlabel('Sun Zenith Angle','FontSize',16,'FontWeight','bold','Color','k')
ylabel('TOA Reflectance','FontSize',16,'FontWeight','bold','Color','k')
title('Before and after BRDF corrected TOA reflectance for CA,Blue,Green & Red','FontSize',20);






hold off
saveas(gcf,'Comparison11'); 
saveas(gcf,['Comparison11','.jpg']); 








figure(2)
for i=5:1:7
 
    subplot(2,2,i-4)
 
    plot(SZA',Refl(:,i),'bo');
    hold on
    plot(SZA',Refl_BRDF_Corrected(:,i),'rd');
    grid on;
    
    
%     
%     v = axis;
% handle=title([Bands{i}])
% set(handle,'Position',[.5 .5]);


  title([Bands{i}],'FontSize',8);
    xlabel('Sun Zenith Angle','FontSize',12,'FontWeight','bold','Color','k')
    ylabel('TOA Reflectance','FontSize',12,'FontWeight','bold','Color','k')
    %legend({'Before','After'},'Location','northwest')
    
    
 ylim([limit1(i) limit2(i)]); 
 
 
 set(gcf,'NextPlot','add');
axes;
h = title('MyTitle');
set(gca,'Visible','off');
set(h,'Visible','on');
 
 
end

% AX=legend('plot 1', 'plot 2','location','northeast','orientation','horizontal');
xlabel('Sun Zenith Angle','FontSize',16,'FontWeight','bold','Color','k')
ylabel('TOA Reflectance','FontSize',16,'FontWeight','bold','Color','k')
title('Before and after BRDF corrected TOA reflectance for NIR,SWIR 1 & SWIR2','FontSize',20);






hold off
saveas(gcf,'Comparison12'); 
saveas(gcf,['Comparison12','.jpg']); 



% 
% for i=1:1:13
%     pppp(i)=str2num(subFolders(i).name)
% end
% 
% 
% T=table(pppp','VariableNames',{'Dates'})
% writetable(T,['Dates' '.csv'])