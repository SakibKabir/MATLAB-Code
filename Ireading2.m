clear all;
r=1100:1200;
c=1200:1300;

rmb=0.00002; % reflectance multiplication band
rab=-0.1; % reflectance add band

%radiance multiplication band
rb1=0.012179;

%radiance add band
rab1= -60.89552;

base='Z:\ImageDrive\OLI.TIRS\L8\P187\R043';
dates = dir(base);
dates([1 2])=[];

for date=1:size(dates,1)
    if exist(fullfile(base,dates(date).name,'LC1'),'dir')
        Band1=dir(fullfile(base,dates(date).name,'LC1','*B1.tif')); % deep blue and violets
        Band2=dir(fullfile(base,dates(date).name,'LC1','*B2.tif')); % visible blue
        Band3=dir(fullfile(base,dates(date).name,'LC1','*B3.tif')); % green
        Band4=dir(fullfile(base,dates(date).name,'LC1','*B4.tif')); % red
        Band5=dir(fullfile(base,dates(date).name,'LC1','*B5.tif')); % NIR
        
        [A1,R1] = geotiffread(fullfile(base,dates(date).name,'LC1',Band1.name));
        dn1(date)=mean2(A1(r,c));
        %ref1(date)=dn1(date)*rmb-rab;
        rad1(date)= dn1(date)*rb1+rab1;
        
        [A2,R2] = geotiffread(fullfile(base,dates(date).name,'LC1',Band2.name));
         dn2(date)=mean2(A2(r,c));
         ref2(date)=dn2(date)*rmb+rab;
         
        % ANew=A1-A2;
         
        [A3,R3] = geotiffread(fullfile(base,dates(date).name,'LC1',Band3.name));
        dn3(date)=mean2(A3(r,c));
        ref3(date)=dn3(date)*rmb+rab;
        
        [A4,R4] = geotiffread(fullfile(base,dates(date).name,'LC1',Band4.name));
        dn4(date)=mean2(A4(r,c));
        ref4(date)=dn4(date)*rmb+rab;
        
        [A5,R5] = geotiffread(fullfile(base,dates(date).name,'LC1',Band5.name));
         dn5(date)=mean2(A5(r,c));
         ref5(date)=dn5(date)*rmb+rab;

    end
    
end

% Newband=ref1-ref2;


subplot(5,1,1);
plot(rad1)

% subplot(5,1,2);
% plot(ref2)
% 
% subplot(5,1,3);
% plot(ref3)
% 
% subplot(5,1,4);
% plot(ref4)
% 
% subplot(5,1,5);
% plot(ref5)
