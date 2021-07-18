

%% Dove R (105e) vs Landsat 7
%%% Landsat L7 Data
leg_location ={'Imperial Valley CA-12/24/18-Scene 51','Imperial Valley CA-12/24/18-Scene 56', 'Imperial Valley CA-12/24/18-Scene 61',...
               'Imperial Valley CA-12/24/18-Scene 66','Imperial Valley CA-12/24/18-Scene 71','Imperial Valley CA-12/24/18-Scene 76'};

load('all_SNR_DR_105e.mat')
load('all_SD_DR_105e.mat')
load('all_mean_DR_105e.mat')

all_mean_rad_L7 = [all_mean_DR_105e.Radiance.a_IV_P038R037_181224_51.L7.ROI; all_mean_DR_105e.Radiance.a_IV_P038R037_181224_56.L7.ROI;...
                   all_mean_DR_105e.Radiance.a_IV_P038R038_181224_61.L7.ROI; all_mean_DR_105e.Radiance.a_IV_P038R038_181224_66.L7.ROI;...
                   all_mean_DR_105e.Radiance.a_IV_P038R038_181224_71.L7.ROI; all_mean_DR_105e.Radiance.a_IV_P038R038_181224_76.L7.ROI]; 
        
all_mean_ref_L7 = [all_mean_DR_105e.Reflectance.a_IV_P038R037_181224_51.L7.ROI; all_mean_DR_105e.Reflectance.a_IV_P038R037_181224_56.L7.ROI;...
                   all_mean_DR_105e.Reflectance.a_IV_P038R038_181224_61.L7.ROI; all_mean_DR_105e.Reflectance.a_IV_P038R038_181224_66.L7.ROI;
                   all_mean_DR_105e.Reflectance.a_IV_P038R038_181224_71.L7.ROI; all_mean_DR_105e.Reflectance.a_IV_P038R038_181224_76.L7.ROI];

%%% Dove-R Data
all_mean_rad_DR = [all_mean_DR_105e.Radiance.a_IV_P038R037_181224_51.D.ROI; all_mean_DR_105e.Radiance.a_IV_P038R037_181224_56.D.ROI;...
                   all_mean_DR_105e.Radiance.a_IV_P038R038_181224_61.D.ROI; all_mean_DR_105e.Radiance.a_IV_P038R038_181224_66.D.ROI;...
                   all_mean_DR_105e.Radiance.a_IV_P038R038_181224_71.D.ROI; all_mean_DR_105e.Radiance.a_IV_P038R038_181224_76.D.ROI];  
             
all_mean_ref_DR = [all_mean_DR_105e.Reflectance.a_IV_P038R037_181224_51.D.ROI; all_mean_DR_105e.Reflectance.a_IV_P038R037_181224_56.D.ROI;...
                   all_mean_DR_105e.Reflectance.a_IV_P038R038_181224_61.D.ROI; all_mean_DR_105e.Reflectance.a_IV_P038R038_181224_66.D.ROI;...
                   all_mean_DR_105e.Reflectance.a_IV_P038R038_181224_71.D.ROI; all_mean_DR_105e.Reflectance.a_IV_P038R038_181224_76.D.ROI];
               
all_SD_rad_DR = [all_SD_DR_105e.Radiance.a_IV_P038R037_181224_51.D.ROI; all_SD_DR_105e.Radiance.a_IV_P038R037_181224_56.D.ROI;...
                 all_SD_DR_105e.Radiance.a_IV_P038R038_181224_61.D.ROI; all_SD_DR_105e.Radiance.a_IV_P038R038_181224_66.D.ROI;...
                 all_SD_DR_105e.Radiance.a_IV_P038R038_181224_71.D.ROI; all_SD_DR_105e.Radiance.a_IV_P038R038_181224_76.D.ROI];            
  
all_SD_ref_DR = [all_SD_DR_105e.Reflectance.a_IV_P038R037_181224_51.D.ROI; all_SD_DR_105e.Reflectance.a_IV_P038R037_181224_56.D.ROI;...
                 all_SD_DR_105e.Reflectance.a_IV_P038R038_181224_61.D.ROI; all_SD_DR_105e.Reflectance.a_IV_P038R038_181224_66.D.ROI;...
                 all_SD_DR_105e.Reflectance.a_IV_P038R038_181224_71.D.ROI; all_SD_DR_105e.Reflectance.a_IV_P038R038_181224_76.D.ROI];

%% L7 vs 1068 Plot
leg_location ={'Imperial Valley CA-01/25/19-Scene 20','Imperial Valley CA-01/25/19-Scene 29', 'Imperial Valley CA-01/25/19-Scene 11',...
               'Imperial Valley CA-01/25/19-Scene 02','Imperial Valley CA-01/25/19-Scene 93', 'Imperial Valley CA-01/25/19-Scene 85'};

load('all_SNR_DR_1068.mat')
load('all_SD_DR_1068.mat')
load('all_mean_DR_1068.mat')

%%% Landsat L7 Data
all_mean_rad_L8 = [all_mean_DR_1068.Radiance.a_IV_P038R037_190125_20.L7.ROI; all_mean_DR_1068.Radiance.a_IV_P038R037_190125_29.L7.ROI;...
                   all_mean_DR_1068.Radiance.a_IV_P038R038_190125_11.L7.ROI; all_mean_DR_1068.Radiance.a_IV_P038R038_190125_02.L7.ROI;...
                   all_mean_DR_1068.Radiance.a_IV_P038R038_190125_93.L7.ROI; all_mean_DR_1068.Radiance.a_IV_P038R038_190125_85.L7.ROI]; 
        
all_mean_ref_L8 = [all_mean_DR_1068.Reflectance.a_IV_P038R037_190125_20.L7.ROI; all_mean_DR_1068.Reflectance.a_IV_P038R037_190125_29.L7.ROI;...
                   all_mean_DR_1068.Reflectance.a_IV_P038R038_190125_11.L7.ROI; all_mean_DR_1068.Reflectance.a_IV_P038R038_190125_02.L7.ROI;
                   all_mean_DR_1068.Reflectance.a_IV_P038R038_190125_93.L7.ROI; all_mean_DR_1068.Reflectance.a_IV_P038R038_190125_85.L7.ROI];

 all_SNR_rad_L8 = [all_SNR_DR_1068.Radiance.a_IV_P038R037_190125_20.L7.ROI; all_SNR_DR_1068.Radiance.a_IV_P038R037_190125_29.L7.ROI;...
                   all_SNR_DR_1068.Radiance.a_IV_P038R038_190125_11.L7.ROI; all_SNR_DR_1068.Radiance.a_IV_P038R038_190125_02.L7.ROI;...
                   all_SNR_DR_1068.Radiance.a_IV_P038R038_190125_93.L7.ROI; all_SNR_DR_1068.Radiance.a_IV_P038R038_190125_85.L7.ROI];               
               
%%% Dove-R Data
all_mean_rad_DR = [all_mean_DR_1068.Radiance.a_IV_P038R037_190125_20.D.ROI; all_mean_DR_1068.Radiance.a_IV_P038R037_190125_29.D.ROI;...
                   all_mean_DR_1068.Radiance.a_IV_P038R038_190125_11.D.ROI; all_mean_DR_1068.Radiance.a_IV_P038R038_190125_02.D.ROI;...
                   all_mean_DR_1068.Radiance.a_IV_P038R038_190125_93.D.ROI; all_mean_DR_1068.Radiance.a_IV_P038R038_190125_85.D.ROI];  
             
all_mean_ref_DR = [all_mean_DR_1068.Reflectance.a_IV_P038R037_190125_20.D.ROI; all_mean_DR_1068.Reflectance.a_IV_P038R037_190125_29.D.ROI;
                   all_mean_DR_1068.Reflectance.a_IV_P038R038_190125_11.D.ROI; all_mean_DR_1068.Reflectance.a_IV_P038R038_190125_02.D.ROI;...
                   all_mean_DR_1068.Reflectance.a_IV_P038R038_190125_93.D.ROI; all_mean_DR_1068.Reflectance.a_IV_P038R038_190125_85.D.ROI];
               
all_SD_rad_DR = [all_SD_DR_1068.Radiance.a_IV_P038R037_190125_20.D.ROI; all_SD_DR_1068.Radiance.a_IV_P038R037_190125_29.D.ROI;...
                 all_SD_DR_1068.Radiance.a_IV_P038R038_190125_11.D.ROI; all_SD_DR_1068.Radiance.a_IV_P038R038_190125_02.D.ROI;...
                 all_SD_DR_1068.Radiance.a_IV_P038R038_190125_93.D.ROI; all_SD_DR_1068.Radiance.a_IV_P038R038_190125_85.D.ROI];            
  
all_SD_ref_DR = [all_SD_DR_1068.Reflectance.a_IV_P038R037_190125_20.D.ROI; all_SD_DR_1068.Reflectance.a_IV_P038R037_190125_29.D.ROI;...
                 all_SD_DR_1068.Reflectance.a_IV_P038R038_190125_11.D.ROI; all_SD_DR_1068.Reflectance.a_IV_P038R038_190125_02.D.ROI;
                 all_SD_DR_1068.Reflectance.a_IV_P038R038_190125_93.D.ROI; all_SD_DR_1068.Reflectance.a_IV_P038R038_190125_85.D.ROI];

%%
%% Dove-R 1058 and Landsat 8
leg_location ={'Imperial Valley CA-02/18/19-Scene 16', 'Imperial Valley CA-02/18/19-Scene 22', 'Imperial Valley CA-02/18/19-Scene 29',...
               'Imperial Valley CA-02/18/19-Scene 36', 'Imperial Valley CA-02/18/19-Scene 42', 'Imperial Valley CA-02/18/19-Scene 49',...
               'Imperial Valley CA-02/18/19-Scene 55'};
           
load('all_mean_DR_1058.mat')
load('all_SD_DR_1058.mat')
load('all_SNR_DR_1058.mat')

%%% Landsat L8 Data
all_mean_rad_L8 = [all_mean_DR_1058.Radiance.a_IV_P038R037_190218_16.L8.ROI; all_mean_DR_1058.Radiance.a_IV_P038R037_190218_22.L8.ROI;...
                   all_mean_DR_1058.Radiance.a_IV_P038R037_190218_29.L8.ROI; all_mean_DR_1058.Radiance.a_IV_P038R038_190218_36.L8.ROI;...
                   all_mean_DR_1058.Radiance.a_IV_P038R038_190218_42.L8.ROI; all_mean_DR_1058.Radiance.a_IV_P038R038_190218_49.L8.ROI;
                   all_mean_DR_1058.Radiance.a_IV_P038R038_190218_55.L8.ROI]; 
        
all_mean_ref_L8 = [all_mean_DR_1058.Reflectance.a_IV_P038R037_190218_16.L8.ROI; all_mean_DR_1058.Reflectance.a_IV_P038R037_190218_22.L8.ROI;...
                   all_mean_DR_1058.Reflectance.a_IV_P038R037_190218_29.L8.ROI; all_mean_DR_1058.Reflectance.a_IV_P038R038_190218_36.L8.ROI;
                   all_mean_DR_1058.Reflectance.a_IV_P038R038_190218_42.L8.ROI; all_mean_DR_1058.Reflectance.a_IV_P038R038_190218_49.L8.ROI;
                   all_mean_DR_1058.Reflectance.a_IV_P038R038_190218_55.L8.ROI];

all_SD_rad_L8 = [all_SD_DR_1058.Radiance.a_IV_P038R037_190218_16.L8.ROI; all_SD_DR_1058.Radiance.a_IV_P038R037_190218_22.L8.ROI;...
                 all_SD_DR_1058.Radiance.a_IV_P038R037_190218_29.L8.ROI; all_SD_DR_1058.Radiance.a_IV_P038R038_190218_36.L8.ROI;
                 all_SD_DR_1058.Radiance.a_IV_P038R038_190218_42.L8.ROI; all_SD_DR_1058.Radiance.a_IV_P038R038_190218_49.L8.ROI;
                 all_SD_DR_1058.Radiance.a_IV_P038R038_190218_55.L8.ROI];
             
all_SD_ref_L8 = [all_SD_DR_1058.Reflectance.a_IV_P038R037_190218_16.L8.ROI; all_SD_DR_1058.Reflectance.a_IV_P038R037_190218_22.L8.ROI;...
                 all_SD_DR_1058.Reflectance.a_IV_P038R037_190218_29.L8.ROI; all_SD_DR_1058.Reflectance.a_IV_P038R038_190218_36.L8.ROI;
                 all_SD_DR_1058.Reflectance.a_IV_P038R038_190218_42.L8.ROI; all_SD_DR_1058.Reflectance.a_IV_P038R038_190218_49.L8.ROI;
                 all_SD_DR_1058.Reflectance.a_IV_P038R038_190218_55.L8.ROI]; 
             
all_SNR_rad_L8 = [all_SNR_DR_1058.Radiance.a_IV_P038R037_190218_16.L8.ROI; all_SNR_DR_1058.Radiance.a_IV_P038R037_190218_22.L8.ROI;...
                 all_SNR_DR_1058.Radiance.a_IV_P038R037_190218_29.L8.ROI; all_SNR_DR_1058.Radiance.a_IV_P038R038_190218_36.L8.ROI;
                 all_SNR_DR_1058.Radiance.a_IV_P038R038_190218_42.L8.ROI; all_SNR_DR_1058.Radiance.a_IV_P038R038_190218_49.L8.ROI;
                 all_SNR_DR_1058.Radiance.a_IV_P038R038_190218_55.L8.ROI];
               
%%% Dove-R Data
all_mean_rad_DR = [all_mean_DR_1058.Radiance.a_IV_P038R037_190218_16.D.ROI; all_mean_DR_1058.Radiance.a_IV_P038R037_190218_22.D.ROI;...
                   all_mean_DR_1058.Radiance.a_IV_P038R037_190218_29.D.ROI; all_mean_DR_1058.Radiance.a_IV_P038R038_190218_36.D.ROI;
                   all_mean_DR_1058.Radiance.a_IV_P038R038_190218_42.D.ROI; all_mean_DR_1058.Radiance.a_IV_P038R038_190218_49.D.ROI;
                   all_mean_DR_1058.Radiance.a_IV_P038R038_190218_55.D.ROI];  
             
all_mean_ref_DR = [all_mean_DR_1058.Reflectance.a_IV_P038R037_190218_16.D.ROI; all_mean_DR_1058.Reflectance.a_IV_P038R037_190218_22.D.ROI;...
                   all_mean_DR_1058.Reflectance.a_IV_P038R037_190218_29.D.ROI; all_mean_DR_1058.Reflectance.a_IV_P038R038_190218_36.D.ROI;
                   all_mean_DR_1058.Reflectance.a_IV_P038R038_190218_42.D.ROI; all_mean_DR_1058.Reflectance.a_IV_P038R038_190218_49.D.ROI;
                   all_mean_DR_1058.Reflectance.a_IV_P038R038_190218_55.D.ROI];
               
all_SD_rad_DR = [all_SD_DR_1058.Radiance.a_IV_P038R037_190218_16.D.ROI; all_SD_DR_1058.Radiance.a_IV_P038R037_190218_22.D.ROI;...
                 all_SD_DR_1058.Radiance.a_IV_P038R037_190218_29.D.ROI; all_SD_DR_1058.Radiance.a_IV_P038R038_190218_36.D.ROI;
                 all_SD_DR_1058.Radiance.a_IV_P038R038_190218_42.D.ROI; all_SD_DR_1058.Radiance.a_IV_P038R038_190218_49.D.ROI;
                 all_SD_DR_1058.Radiance.a_IV_P038R038_190218_55.D.ROI];            

all_SD_ref_DR = [all_SD_DR_1058.Reflectance.a_IV_P038R037_190218_16.D.ROI; all_SD_DR_1058.Reflectance.a_IV_P038R037_190218_22.D.ROI;
                 all_SD_DR_1058.Reflectance.a_IV_P038R037_190218_29.D.ROI; all_SD_DR_1058.Reflectance.a_IV_P038R038_190218_36.D.ROI;
                 all_SD_DR_1058.Reflectance.a_IV_P038R038_190218_42.D.ROI; all_SD_DR_1058.Reflectance.a_IV_P038R038_190218_49.D.ROI;
                 all_SD_DR_1058.Reflectance.a_IV_P038R038_190218_55.D.ROI];

all_SNR_rad_DR =[all_SNR_DR_1058.Radiance.a_IV_P038R037_190218_16.D.ROI; all_SNR_DR_1058.Radiance.a_IV_P038R037_190218_22.D.ROI;...
                 all_SNR_DR_1058.Radiance.a_IV_P038R037_190218_29.D.ROI; all_SNR_DR_1058.Radiance.a_IV_P038R038_190218_36.D.ROI;
                 all_SNR_DR_1058.Radiance.a_IV_P038R038_190218_42.D.ROI; all_SNR_DR_1058.Radiance.a_IV_P038R038_190218_49.D.ROI;
                 all_SNR_DR_1058.Radiance.a_IV_P038R038_190218_55.D.ROI];
             
%%
figure, plot (all_mean_rad_L8(:, 1),all_mean_rad_DR(:,1), '.' )

figure, plot (all_mean_rad_L8(:, 1), all_SNR_rad_L8(:, 1),  '.' )

figure, plot (all_mean_DR_1058.Radiance.a_IV_P038R038_190218_36.D.ROI(:,1), all_SNR_DR_1058.Radiance.a_IV_P038R038_190218_36.D.ROI(:,1), '.')

all_mean_rad_DR
all_SNR_rad_DR

xlim([0 200])

