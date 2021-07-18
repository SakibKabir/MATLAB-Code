close all

figure
image_na=imread('LC08_L1TP_161055_20200407_20200410_01_T1_B1.tif');
imagesc(image_na), colormap(gray);

figure
esfh = image_na(2450:2490, 3780:3835);
imagesc(esfh),colormap(gray)
%%
for line = 1:50
    plot(esfh (line, :),'linestyle','none','MarkerSize',15,'Marker','.'); hold on 
    title('Horizontal edge spread function ESP for 50 scans'); 
    xlabel('index number')
    ylabel('pixels (rows)')
    grid on;  grid minor; 
    set(gca, 'FontSize', 15, 'Fontweight','bold'); set (gcf,'color','white'); 
end 

%%
mean_esfh=mean(esfh);
figure 
plot(mean_esfh,'linestyle','none','MarkerSize',25,'Marker','.','color','m');  
title('Mean of the horizontal edge spread function ESP for 50 scans');
xlabel('Index number') 
ylabel('pixels (rows)')
grid on;  grid minor; 
set(gca, 'FontSize', 15, 'Fontweight','bold'); 
set (gcf,'color','white'); 

LSFH1 = mean_esfh;
LSFH1 = double(LSFH1);
LSFH = sgolayfilt(LSFH1,3,11);
figure
plot(LSFH);

%%
a_downh = mean(LSFH (1,37:56)); 
a_uph = mean(LSFH (1,1:16)); 
a_totalh = abs(a_downh - a_uph ); 
bh= LSFH (27); ch =19; 
dh = a_downh;
fermi_h = a_totalh./(exp(-(LSFH -bh)./ch)+1)+dh;
figure, 
plot(LSFH , 'linestyle','none','MarkerSize',30,'Marker','.','color','m');hold on 
plot(fermi_h, 'MarkerSize',30,'Marker','.','color','k');hold on  
title('Mean of the horizontal edge spread function ESP for 50 scans and the fitted Fermi-function'); 
xlabel('Index number') 
ylabel('pixels (rows)') 
grid on;  

%%
LSF_H = diff(fermi_h); 
figure 
plot(LSF_H, 'linestyle','-','MarkerSize',25,'Marker','.','color','m');  
title('Horizontal line spread function (LSF)');
xlabel('index number')
ylabel('pixels (rows)') 
grid on; 
grid minor; 
set(gca, 'FontSize', 15, 'Fontweight','bold');
set (gcf,'color','white'); 
%%
MTF = abs(fft(LSF_H))./max(abs(fft(LSF_H)));
figure
plot (MTF,'.')


