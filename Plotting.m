%% Ploting for all the bands
close all
band_name = {'BLUE','GREEN' ,'RED' ,'NIR' };
band_colors={'c','g','r','m'};
leg ={'ROI 1', 'ROI 2', 'ROI 3',  'ROI 4', 'ROI 5',  'ROI 6', 'ROI 7'};
leg_location ={'California', 'Wisconsin', 'Indiana',  'Missouri', 'J,SA(P/R-170/078)', 'J,SA(P/R-170/079)', ...
                'J,SA(P/R-171/078)', 'J,SA(P/R-171/079)','CT1,SA(P/R-174/084)', 'CT2,SA(P/R-175/083)',...
               'Auckland,NZ', 'Sydney, Australia', 'NSW,Aus(P/R-090/082)', 'NSW,Aus(P/R-091/082)',...
               'Melbourne, Australia', 'Libya4,Scene-52', 'Libya4,Scene-50'};

all_date = {'2018-12-23', '2018-09-23', '2018-11-07', '2018-11-21', '2018-12-13', '2018-12-13', '2018-12-04', '2018-12-04',...
            '2018-09-20', '2018-09-11', '2018-09-09', '2018-11-20', '2018-11-11', '2018-11-02', '2018-11-16', '2018-10-07', '2018-10-07'};           

location_str = string(leg_location);

%location = {'California', 'Wisconsin', 'Indiana',  'Missouri', 'J-SA'};
marker = {'o', 's', 'd', '+','p',  '^','>', '<', 'v', 'o', 's', 'd', '+','p','x','*','h', '+', '.'};

%%%% Radiance
%%% Mean
    allmean_rad_L8 = [all_mean_com.Radiance.a_california.L8.ROI; all_mean_com.Radiance.b_wisconsin.L8.ROI;...
                      all_mean_com.Radiance.c_indiana.L8.ROI; all_mean_com.Radiance.d_missouri.L8.ROI;...
                      all_mean_com.Radiance.e_J_SA1.L8.ROI; all_mean_com.Radiance.f_J_SA2.L8.ROI;...
                      all_mean_com.Radiance.g_J_SA3.L8.ROI; all_mean_com.Radiance.h_J_SA4.L8.ROI;...
                      all_mean_com.Radiance.i_CT_SA1.L8.ROI; all_mean_com.Radiance.j_CT_SA2.L8.ROI;...
                      all_mean_com.Radiance.k_A_NZ.L8.ROI; all_mean_com.Radiance.l_S_Aus.L8.ROI;...
                      all_mean_com.Radiance.m_NSW_Aus1.L8.ROI; all_mean_com.Radiance.n_NSW_Aus2.L8.ROI;...
                      all_mean_com.Radiance.o_MB_Aus.L8.ROI; all_mean_com.Radiance.p_Libya4_52.L8.ROI;...
                      all_mean_com.Radiance.p_Libya4_50.L8.ROI];
                                   
    allmean_rad_D = [all_mean_com.Radiance.a_california.D.ROI; all_mean_com.Radiance.b_wisconsin.D.ROI;...
                     all_mean_com.Radiance.c_indiana.D.ROI; all_mean_com.Radiance.d_missouri.D.ROI;...
                     all_mean_com.Radiance.e_J_SA1.D.ROI; all_mean_com.Radiance.f_J_SA2.D.ROI;...
                     all_mean_com.Radiance.g_J_SA3.D.ROI; all_mean_com.Radiance.h_J_SA4.D.ROI;...
                     all_mean_com.Radiance.i_CT_SA1.D.ROI; all_mean_com.Radiance.j_CT_SA2.D.ROI;...
                     all_mean_com.Radiance.k_A_NZ.D.ROI; all_mean_com.Radiance.l_S_Aus.D.ROI;...
                     all_mean_com.Radiance.m_NSW_Aus1.D.ROI; all_mean_com.Radiance.n_NSW_Aus2.D.ROI;...
                     all_mean_com.Radiance.o_MB_Aus.D.ROI; all_mean_com.Radiance.p_Libya4_52.D.ROI;...
                     all_mean_com.Radiance.p_Libya4_50.D.ROI];
                 
%    allmean_rad_D = [all_mean_com.Radiance.a_california.D.ROI; all_mean_com.Radiance.b_wisconsin.D.ROI;...
%                      all_mean_com.Radiance.c_indiana.D.ROI; all_mean_com.Radiance.d_missouri.D.ROI;...
%                      all_mean.Radiance.e_J_SA1.D.ROI; all_mean_com.Radiance.f_J_SA2.D.ROI;...
%                      all_mean.Radiance.g_J_SA3.D.ROI; all_mean_com.Radiance.h_J_SA4.D.ROI;...
%                      all_mean_com.Radiance.i_CT_SA1.D.ROI; all_mean_com.Radiance.j_CT_SA2.D.ROI;...
%                      all_mean_com.Radiance.k_A_NZ.D.ROI; all_mean_com.Radiance.l_S_Aus.D.ROI;...
%                      all_mean_com.Radiance.m_NSW_Aus1.D.ROI; all_mean_com.Radiance.n_NSW_Aus2.D.ROI;...
%                      all_mean_com.Radiance.o_MB_Aus.D.ROI];
   
    %%% Reflectance   
    allmean_ref_L8 = [all_mean_com.Reflectance.a_california.L8.ROI; all_mean_com.Reflectance.b_wisconsin.L8.ROI;...
                     all_mean_com.Reflectance.c_indiana.L8.ROI; all_mean_com.Reflectance.d_missouri.L8.ROI;...
                     all_mean_com.Reflectance.e_J_SA1.L8.ROI; all_mean_com.Reflectance.f_J_SA2.L8.ROI;...
                     all_mean_com.Reflectance.g_J_SA3.L8.ROI; all_mean_com.Reflectance.h_J_SA4.L8.ROI;...
                     all_mean_com.Reflectance.i_CT_SA1.L8.ROI; all_mean_com.Reflectance.j_CT_SA2.L8.ROI;...
                     all_mean_com.Reflectance.k_A_NZ.L8.ROI; all_mean_com.Reflectance.l_S_Aus.L8.ROI;...
                     all_mean_com.Reflectance.m_NSW_Aus1.L8.ROI; all_mean_com.Reflectance.n_NSW_Aus2.L8.ROI;...
                     all_mean_com.Reflectance.o_MB_Aus.L8.ROI; all_mean_com.Reflectance.p_Libya4_52.L8.ROI;...
                     all_mean_com.Reflectance.p_Libya4_50.L8.ROI];
   
    allmean_ref_D = [all_mean_com.Reflectance.a_california.D.ROI; all_mean_com.Reflectance.b_wisconsin.D.ROI;...
                     all_mean_com.Reflectance.c_indiana.D.ROI; all_mean_com.Reflectance.d_missouri.D.ROI;...
                     all_mean_com.Reflectance.e_J_SA1.D.ROI; all_mean_com.Reflectance.f_J_SA2.D.ROI;...
                     all_mean_com.Reflectance.g_J_SA3.D.ROI; all_mean_com.Reflectance.h_J_SA4.D.ROI;...
                     all_mean_com.Reflectance.i_CT_SA1.D.ROI; all_mean_com.Reflectance.j_CT_SA2.D.ROI;...
                     all_mean_com.Reflectance.k_A_NZ.D.ROI; all_mean_com.Reflectance.l_S_Aus.D.ROI;....
                     all_mean_com.Reflectance.m_NSW_Aus1.D.ROI; all_mean_com.Reflectance.n_NSW_Aus2.D.ROI;...
                     all_mean_com.Reflectance.o_MB_Aus.D.ROI; all_mean_com.Reflectance.p_Libya4_52.D.ROI;
                     all_mean_com.Reflectance.p_Libya4_50.D.ROI];
   
  %%% Relative Standard Deviation
   % Radiance
    all_RSD_rad_L8 = [all_RSD.Radiance.a_california.L8.ROI; all_RSD.Radiance.b_wisconsin.L8.ROI;...
                     all_RSD.Radiance.c_indiana.L8.ROI; all_RSD.Radiance.d_missouri.L8.ROI;...
                     all_RSD.Radiance.e_J_SA1.L8.ROI; all_RSD.Radiance.f_J_SA2.L8.ROI;...
                     all_RSD.Radiance.g_J_SA3.L8.ROI; all_RSD.Radiance.h_J_SA4.L8.ROI;...
                     all_RSD.Radiance.i_CT_SA1.L8.ROI; all_RSD.Radiance.j_CT_SA2.L8.ROI;...
                     all_RSD.Radiance.k_A_NZ.L8.ROI; all_RSD.Radiance.l_S_Aus.L8.ROI;...
                     all_RSD.Radiance.m_NSW_Aus1.L8.ROI; all_RSD.Radiance.n_NSW_Aus2.L8.ROI;...
                     all_RSD.Radiance.o_MB_Aus.L8.ROI; all_RSD.Radiance.p_Libya4_52.L8.ROI;...
                     all_RSD.Radiance.p_Libya4_50.L8.ROI];
   
     all_RSD_rad_D =[all_RSD.Radiance.a_california.D.ROI; all_RSD.Radiance.b_wisconsin.D.ROI;...
                    all_RSD.Radiance.c_indiana.D.ROI; all_RSD.Radiance.d_missouri.D.ROI;...
                    all_RSD.Radiance.e_J_SA1.D.ROI; all_RSD.Radiance.f_J_SA2.D.ROI;...
                    all_RSD.Radiance.g_J_SA3.D.ROI; all_RSD.Radiance.h_J_SA4.D.ROI;...
                    all_RSD.Radiance.i_CT_SA1.D.ROI; all_RSD.Radiance.j_CT_SA2.D.ROI;...
                    all_RSD.Radiance.k_A_NZ.D.ROI; all_RSD.Radiance.l_S_Aus.D.ROI;...
                    all_RSD.Radiance.m_NSW_Aus1.D.ROI; all_RSD.Radiance.n_NSW_Aus2.D.ROI;...
                    all_RSD.Radiance.o_MB_Aus.D.ROI; all_RSD.Radiance.p_Libya4_52.D.ROI;...
                    all_RSD.Radiance.p_Libya4_50.D.ROI];
   
% Reflectance
    all_RSD_ref_L8 =[all_RSD.Reflectance.a_california.L8.ROI; all_RSD.Reflectance.b_wisconsin.L8.ROI;...
                     all_RSD.Reflectance.c_indiana.L8.ROI; all_RSD.Reflectance.d_missouri.L8.ROI;...
                     all_RSD.Reflectance.e_J_SA1.L8.ROI; all_RSD.Reflectance.f_J_SA2.L8.ROI;...
                     all_RSD.Reflectance.g_J_SA3.L8.ROI; all_RSD.Reflectance.h_J_SA4.L8.ROI;...
                     all_RSD.Reflectance.i_CT_SA1.L8.ROI; all_RSD.Reflectance.j_CT_SA2.L8.ROI;...
                     all_RSD.Reflectance.k_A_NZ.L8.ROI; all_RSD.Reflectance.l_S_Aus.L8.ROI;...
                     all_RSD.Reflectance.m_NSW_Aus1.L8.ROI; all_RSD.Reflectance.n_NSW_Aus2.L8.ROI;...
                     all_RSD.Reflectance.o_MB_Aus.L8.ROI; all_RSD.Reflectance.p_Libya4_52.L8.ROI;...
                     all_RSD.Reflectance.p_Libya4_50.L8.ROI];
   
     all_RSD_ref_D =[all_RSD.Reflectance.a_california.D.ROI; all_RSD.Reflectance.b_wisconsin.D.ROI;...
                   all_RSD.Reflectance.c_indiana.D.ROI; all_RSD.Reflectance.d_missouri.D.ROI;...
                   all_RSD.Reflectance.e_J_SA1.D.ROI; all_RSD.Reflectance.f_J_SA2.D.ROI;...
                   all_RSD.Reflectance.g_J_SA3.D.ROI; all_RSD.Reflectance.h_J_SA4.D.ROI;...
                   all_RSD.Reflectance.i_CT_SA1.D.ROI; all_RSD.Reflectance.j_CT_SA2.D.ROI;...
                   all_RSD.Reflectance.k_A_NZ.D.ROI; all_RSD.Reflectance.l_S_Aus.D.ROI;...
                   all_RSD.Reflectance.m_NSW_Aus1.D.ROI; all_RSD.Reflectance.n_NSW_Aus2.D.ROI;...
                   all_RSD.Reflectance.o_MB_Aus.D.ROI; all_RSD.Reflectance.p_Libya4_52.D.ROI;
                   all_RSD.Reflectance.p_Libya4_50.D.ROI];

%%% Standard Deviation
    % Radiance
    all_SD_rad_L8 = [all_SD.Radiance.a_california.L8.ROI; all_SD.Radiance.b_wisconsin.L8.ROI;...
                      all_SD.Radiance.c_indiana.L8.ROI; all_SD.Radiance.d_missouri.L8.ROI;...
                      all_SD.Radiance.e_J_SA1.L8.ROI; all_SD.Radiance.f_J_SA2.L8.ROI;...
                      all_SD.Radiance.g_J_SA3.L8.ROI; all_SD.Radiance.h_J_SA4.L8.ROI;...
                      all_SD.Radiance.i_CT_SA1.L8.ROI; all_SD.Radiance.j_CT_SA2.L8.ROI;...
                      all_SD.Radiance.k_A_NZ.L8.ROI; all_SD.Radiance.l_S_Aus.L8.ROI;...
                      all_SD.Radiance.m_NSW_Aus1.L8.ROI; all_SD.Radiance.n_NSW_Aus2.L8.ROI;...
                      all_SD.Radiance.o_MB_Aus.L8.ROI; all_SD.Radiance.p_Libya4_52.L8.ROI;...
                      all_SD.Radiance.p_Libya4_52.L8.ROI];

    all_SD_rad_D = [all_SD.Radiance.a_california.D.ROI; all_SD.Radiance.b_wisconsin.D.ROI;...
                     all_SD.Radiance.c_indiana.D.ROI; all_SD.Radiance.d_missouri.D.ROI;...
                     all_SD.Radiance.e_J_SA1.D.ROI; all_SD.Radiance.f_J_SA2.D.ROI;...
                     all_SD.Radiance.g_J_SA3.D.ROI; all_SD.Radiance.h_J_SA4.D.ROI;...
                     all_SD.Radiance.i_CT_SA1.D.ROI; all_SD.Radiance.j_CT_SA2.D.ROI;...
                     all_SD.Radiance.k_A_NZ.D.ROI; all_SD.Radiance.l_S_Aus.D.ROI;...
                     all_SD.Radiance.m_NSW_Aus1.D.ROI; all_SD.Radiance.n_NSW_Aus2.D.ROI;...
                     all_SD.Radiance.o_MB_Aus.D.ROI; all_SD.Radiance.p_Libya4_52.D.ROI;...
                     all_SD.Radiance.p_Libya4_50.D.ROI;];
                 
   % Reflectance
    all_SD_ref_L8 = [all_SD.Reflectance.a_california.L8.ROI; all_SD.Reflectance.b_wisconsin.L8.ROI;...
                      all_SD.Reflectance.c_indiana.L8.ROI; all_SD.Reflectance.d_missouri.L8.ROI;...
                      all_SD.Reflectance.e_J_SA1.L8.ROI; all_SD.Reflectance.f_J_SA2.L8.ROI;...
                      all_SD.Reflectance.g_J_SA3.L8.ROI; all_SD.Reflectance.h_J_SA4.L8.ROI;...
                      all_SD.Reflectance.i_CT_SA1.L8.ROI; all_SD.Reflectance.j_CT_SA2.L8.ROI;...
                      all_SD.Reflectance.k_A_NZ.L8.ROI; all_SD.Reflectance.l_S_Aus.L8.ROI;...
                      all_SD.Reflectance.m_NSW_Aus1.L8.ROI; all_SD.Reflectance.n_NSW_Aus2.L8.ROI;...
                      all_SD.Reflectance.o_MB_Aus.L8.ROI; all_SD.Reflectance.p_Libya4_52.L8.ROI;...
                      all_SD.Reflectance.p_Libya4_50.L8.ROI;];

    all_SD_ref_D = [all_SD.Reflectance.a_california.D.ROI; all_SD.Reflectance.b_wisconsin.D.ROI;...
                     all_SD.Reflectance.c_indiana.D.ROI; all_SD.Reflectance.d_missouri.D.ROI;...
                     all_SD.Reflectance.e_J_SA1.D.ROI; all_SD.Reflectance.f_J_SA2.D.ROI;...
                     all_SD.Reflectance.g_J_SA3.D.ROI; all_SD.Reflectance.h_J_SA4.D.ROI;...
                     all_SD.Reflectance.i_CT_SA1.D.ROI; all_SD.Reflectance.j_CT_SA2.D.ROI;...
                     all_SD.Reflectance.k_A_NZ.D.ROI; all_SD.Reflectance.l_S_Aus.D.ROI;...
                     all_SD.Reflectance.m_NSW_Aus1.D.ROI; all_SD.Reflectance.n_NSW_Aus2.D.ROI;...
                     all_SD.Reflectance.o_MB_Aus.D.ROI; all_SD.Reflectance.p_Libya4_52.D.ROI;
                     all_SD.Reflectance.p_Libya4_50.D.ROI];
    %%
    close all
 for band=1:4
   if band==1    
        % band =1
        %%% Mean
        L8_band2 = allmean_rad_L8(:,band); D_band1 = allmean_rad_D(:,band); %Radiance
        L8_band2_ref = allmean_ref_L8(:,band); D_band1_ref = allmean_ref_D(:,band); %Reflectance
        
        %%% Standard Deviation
        L8_band2_SD = all_SD_rad_L8(:,band); D_band1_SD = all_SD_rad_D(:,band); %Radiance
        L8_band2_ref_SD = all_SD_ref_L8(:,band); D_band1_ref_SD = all_SD_ref_D(:,band); %Reflectance
        
        %%% Relative Standard Deviation
        L8_band2_RSD = all_RSD_rad_L8(:,band); D_band1_RSD = all_RSD_rad_D(:,band); %Radiance
        L8_band2_ref_RSD = all_RSD_ref_L8(:,band); D_band1_ref_RSD = all_RSD_ref_D(:,band); %Reflectance
       
        %%% Signal to Noise Ratio
        SNR_band1_D = D_band1./D_band1_SD; 
        SNR_band2_L8 = L8_band2./L8_band2_SD;
        
        % Reflectance
        SNR_band1_D_ref = D_band1_ref./D_band1_ref_SD;
        SNR_band2_L8_ref = L8_band2_ref./L8_band2_ref_SD; 
        
        sz=50; c=1; cal=4; 
        for l=1:5
            figure(1)
            %subplot(2,2,1)
            % Radiance
            scatter(D_band1(c:cal), SNR_band1_D(c:cal), sz, marker{l},'LineWidth',1)
            %xlim([0 350]); ylim([0 150]);

            % Reflectance
            %scatter(D_band1_ref(c:cal), SNR_band1_D_ref(c:cal), sz, marker{l},'LineWidth',1)
            %scatter(L8_band2_ref(c:cal), SNR_band2_L8_ref(c:cal), sz, marker{l},'LineWidth',1)
            %xlim([0 0.5]); ylim([0 120]);
            c=cal+1; cal=cal+4;  
            hold on
        end 
        
        legend(leg_location,'FontSize',16, 'Location','southeast');
        
    elseif band ==2
        % band=2
        L8_band3 = allmean_rad_L8(:,band); D_band2 = allmean_rad_D(:,band);
        L8_band3_ref = allmean_ref_L8(:,band); D_band2_ref = allmean_ref_D(:,band);
        
        % Standard Deviation
        L8_band3_SD = all_SD_rad_L8(:,band); D_band2_SD = all_SD_rad_D(:,band); %Radiance
        L8_band3_ref_SD = all_SD_ref_L8(:,band); D_band2_ref_SD = all_SD_ref_D(:,band); %Reflectance
        
        % Relative Standard Deviation
        L8_band3_RSD = all_RSD_rad_L8(:,band); D_band2_RSD = all_RSD_rad_D(:,band); %Radiance
        L8_band3_ref_RSD = all_RSD_ref_L8(:,band); D_band2_ref_RSD = all_RSD_ref_D(:,band); %Reflectance
        
        %%% Signal to Noise Ratio
        SNR_band2_D = D_band2./D_band2_SD; SNR_band3_L8 = L8_band3./L8_band3_SD;
        
        % Reflectance
        SNR_band2_D_ref = D_band2_ref./D_band2_ref_SD; SNR_band3_L8_ref = L8_band3_ref./L8_band3_ref_SD;
        
        hold on; 
        sz=50; c=1; cal=4; 
        for l=1:17
            figure(2)
            %subplot(2,2,2)
            scatter(D_band2(c:cal), SNR_band2_D(c:cal), sz, marker{l},'LineWidth',1)
            xlim([0 350]); ylim([0 150]);
            
            % Reflectance
           % scatter(D_band2_ref(c:cal), SNR_band2_D_ref(c:cal), sz, marker{l},'LineWidth',1)
%             scatter(L8_band3_ref(c:cal), SNR_band3_L8_ref(c:cal), sz, marker{l},'LineWidth',1)
%             xlim([0 0.6]); ylim([0 100]);
%             
            c=cal+1; cal=cal+4;  
            hold on
        end 
        legend(leg_location,'FontSize',16, 'Location','southeast');
        
     elseif band == 3
         % band=2
        L8_band4 = allmean_rad_L8(:,band); D_band3 = allmean_rad_D(:,band);
        L8_band4_ref = allmean_ref_L8(:,band); D_band3_ref = allmean_ref_D(:,band);
        
        % Standard Deviation
        L8_band4_SD = all_SD_rad_L8(:,band); D_band3_SD = all_SD_rad_D(:,band); %Radiance
        L8_band4_ref_SD = all_SD_ref_L8(:,band); D_band3_ref_SD = all_SD_ref_D(:,band); %Reflectance
        
        % Relative Standard Deviation
        L8_band4_RSD = all_RSD_rad_L8(:,band); D_band3_RSD = all_RSD_rad_D(:,band); %Radiance
        L8_band4_ref_RSD = all_RSD_ref_L8(:,band); D_band3_ref_RSD = all_RSD_ref_D(:,band); %Reflectance
        
        %%% Signal to Noise Ratio
        SNR_band3_D = D_band3./D_band3_SD; SNR_band4_L8 = L8_band4./L8_band4_SD;
        
        % Reflectance
        SNR_band3_D_ref = D_band3_ref./D_band3_ref_SD; SNR_band4_L8_ref = L8_band4_ref./L8_band4_ref_SD;
        
        sz=50; c=1; cal=4; 
        for l=1:17
            figure(3)
            %subplot(2,2,3)
            % Radiance
            scatter(D_band3(c:cal), SNR_band3_D(c:cal) , sz, marker{l},'LineWidth',1)
            xlim([0 350]); ylim([0 150]);
            
            % Reflectance
            %scatter(D_band3_ref(c:cal), SNR_band3_D_ref(c:cal), sz, marker{l},'LineWidth',1)
%             scatter(L8_band4_ref(c:cal), SNR_band4_L8_ref(c:cal), sz, marker{l},'LineWidth',1)  
%             xlim([0 0.6]); ylim([0 100]);
            
            c=cal+1; cal=cal+4;  
            hold on
        end 
        legend(leg_location,'FontSize',16, 'Location','southeast');
        
    elseif band == 4
        % band=2
        L8_band5 = allmean_rad_L8(:,band); D_band4 = allmean_rad_D(:,band);
        L8_band5_ref = allmean_ref_L8(:,band); D_band4_ref = allmean_ref_D(:,band);
        
        % Standard Deviation
        L8_band5_SD = all_SD_rad_L8(:,band); D_band4_SD = all_SD_rad_D(:,band); %Radiance
        L8_band5_ref_SD = all_SD_ref_L8(:,band); D_band4_ref_SD = all_SD_ref_D(:,band); %Reflectance
        
        % Relative Standard Deviation
        L8_band5_RSD = all_RSD_rad_L8(:,band); D_band4_RSD = all_RSD_rad_D(:,band); %Radiance
        L8_band5_ref_RSD = all_RSD_ref_L8(:,band); D_band4_ref_RSD = all_RSD_ref_D(:,band); %Reflectance
        
        %%% Signal to Noise Ratio
        SNR_band4_D = D_band4./D_band4_SD; SNR_band4_L8 = L8_band5./L8_band5_SD;
        
        % Reflectance
        SNR_band4_D_ref = D_band4_ref./D_band4_ref_SD; SNR_band5_L8_ref = L8_band5_ref./L8_band5_ref_SD;
        
        sz=50; c=1; cal=4; 
        for l=1:17
            figure(4)
            %subplot(2,2,4)
            % Radiance
            scatter(D_band4(c:cal), SNR_band4_D(c:cal) , sz, marker{l},'LineWidth',1)
            xlim([0 350]); ylim([0 150]);
            
            % Reflectance
            %scatter(D_band4_ref(c:cal), SNR_band4_D_ref(c:cal), sz, marker{l},'LineWidth',1)
%             scatter(L8_band5_ref(c:cal), SNR_band5_L8_ref(c:cal), sz, marker{l},'LineWidth',1)
%             xlim([0 1]); ylim([0 120]);
            
            c=cal+1; cal=cal+4;  
            hold on
        end
        legend(leg_location,'FontSize',10, 'Location','southeast');
        
   end
   
   
    %%% Dove   
    % Radiance
    title(strcat('SNR vs. Mean TOA Radiance', {', '}, strcat(band_name{band}, ' Band')));
    xlabel('Mean TOA Radiance 1047 ROI (W/Sr/m^2/{\mum})')
    ylabel('Signal-to-Noise Ratio')
    
    % Reflectance
%     title(strcat('SNR vs. Mean TOA Reflectance', {', '}, strcat(band_name{band}, ' Band')));
%     xlabel('Mean TOA Reflectance of 1047 ROI')
%     ylabel('Signal-to-Noise Ratio')
%     
    %%% Landsat 8
    % Radiance
%     title(strcat('SNR vs. Mean TOA Radiance', {', '}, strcat(band_name{band}, ' Band')));
%     xlabel('Mean TOA Radiance L8 ROI (W/Sr/m^2/{\mum})')
%     ylabel('Signal-to-Noise Ratio')
    
    % Reflectance
%     title(strcat('SNR vs. Mean TOA Reflectance', {', '}, strcat(band_name{band}, ' Band')));
%     xlabel('Mean TOA Reflectance of L8 ROI')
%     ylabel('Signal-to-Noise Ratio')
%     
    hold on
    grid on
    grid minor
    ax  = gca;
    ax.FontSize = 18;
    ax.GridColor = 'k';
%   ax.GridAlpha = 0.8;
    ax.MinorGridColor = 'k';
    %ax.FontName = 'Times New Roman';
 end     
        %%             
 for band=1:4 %Dove band 
   t= 140;
   if band==1
       
        % band =1
        %%% Mean
        L8_band2 = allmean_rad_L8(:,band); D_band1 = allmean_rad_D(:,band); %Radiance
        L8_band2_ref = allmean_ref_L8(:,band); D_band1_ref = allmean_ref_D(:,band); %Reflectance
        
        %%% Standard Deviation
        L8_band2_SD = all_SD_rad_L8(:,band); D_band1_SD = all_SD_rad_D(:,band); %Radiance
        L8_band2_ref_SD = all_SD_ref_L8(:,band); D_band1_ref_SD = all_SD_ref_D(:,band); %Reflectance
        
        %%% Relative Standard Deviation
        L8_band2_RSD = all_RSD_rad_L8(:,band); D_band1_RSD = all_RSD_rad_D(:,band); %Radiance
        L8_band2_ref_RSD = all_RSD_ref_L8(:,band); D_band1_ref_RSD = all_RSD_ref_D(:,band); %Reflectance
       
        %%% Signal to Noise Ratio
        SNR_band1_D = D_band1./D_band1_SD; 
        SNR_band2_L8 = L8_band2./L8_band2_SD;
        
        % Reflectance
        SNR_band1_D_ref = D_band1_ref./D_band1_ref_SD;
        SNR_band2_L8_ref = L8_band2_ref./L8_band2_ref_SD;
        
        sz=50; c=1; cal=4; 
        
        for l=1:15
            % Radiance
%             p_slope = polyfit(L8_band2(c:cal),D_band1(c:cal),1);
%             p_slope = round(p_slope, 5);
%             fitline_ind = fit(L8_band2(c:cal),D_band1(c:cal), 'poly1');
            
            % Reflectance
            p_slope = polyfit(L8_band2_ref(c:cal), D_band1_ref(c:cal), 1);
            p_slope = round(p_slope, 5);
            fitline_ind = fit(L8_band2_ref(c:cal), D_band1_ref(c:cal), 'poly1');
            
            % Display evaluated equation y = m*x + b
            equation = [num2str(p_slope(1)) '*x + ' num2str(p_slope(2))];
            slope_band1(l)= p_slope(1); bias_band1(l)=p_slope(2);
            
            date = all_date(l);
            [doy,fraction] = date2doy(datenum(char(date)));
            DoY(l)=doy+fraction; %day of year

            G_and_B = coeffvalues(fitline_ind);
            Gain_band1(l) = G_and_B(1,1);
            Bias_band1(l) = G_and_B(1,2);
            ci = confint(fitline_ind);
            ci = round(ci, 5);
            Gain_ci_band1= [ci(1,1) ci(2,1)]; 
            Bias_ci_band1= [ci(1,2) ci(2,2)];
            All_Gain.Gain_ci_band1(l,:) = Gain_ci_band1;
            All_Bias.Bias_ci_band1(l,:) = Bias_ci_band1;
            
            % Evaluate fit equation using polyval
%           f = polyval(p_slope, L8_band2(c:cal));
            
            figure(1);
%             f = fit(L8_band2,D_band1,'poly1');
%             plot(f, L8_band2, D_band1)

%             scatter(L8_band2(c:cal), D_band1(c:cal), sz, marker{l})
%             xlim([0 250]); ylim([0 250])
%             scatter(L8_band2_ref(c:cal), D_band1_ref(c:cal), sz, marker{l})
%             xlim([0 0.6]); ylim([0 0.6]);
            
            % Radiance
%             error = D_band1_SD(c:cal);
%             errorbar(L8_band2(c:cal), D_band1(c:cal), error, marker{l},'LineWidth',0.25)
%             xlim([0 250]); ylim([0 250]);
            
%             % Reflectance
            error = D_band1_ref_SD(c:cal);
            errorbar(L8_band2_ref(c:cal), D_band1_ref(c:cal), error, marker{l},'LineWidth',0.25)
            xlim([0 0.6]); ylim([0 0.6]);
            
            hold on
            M = location_str(l);
            c=cal+1;
            cal=cal+4;
%             tx= strcat('y= ', '  ', equation, ' (', char(M), ')');
%             text(10, t, tx, 'FontSize', 12);
%             t= t + 7;
%             fprintf('%s\n', tx)
        end
     
        legend(leg_location,'FontSize',16, 'Location','southeast');
        hold on
       
        %%% Radiance
%         [p_slope, S] = polyfit(L8_band2, D_band1,1);
%         p_slope = round(p_slope, 5);
%         f = polyval(p_slope, L8_band2); plot(L8_band2, f, '-r'); hold on
%         fitline = fit(L8_band2, D_band1, 'poly1');
        
        %%% Reflectance
        [p_slope, S] = polyfit(L8_band2_ref, D_band1_ref,1);
        p_slope = round(p_slope, 5);
        fitline = fit(L8_band2_ref, D_band1_ref, 'poly1');
        f = polyval(p_slope, L8_band2_ref); plot(L8_band2_ref, f); hold on
          
        %%% Confidence bound
        G_and_B = coeffvalues(fitline);
        Gain_band1 = G_and_B(1,1);
        Bias_band1 = G_and_B(1,2);
        ci = confint(fitline);
        ci = round(ci, 5);
        Gain_lower_bound = ci(1,1); Gain_upper_bound = ci(2,1);
        Bias_lower_bound = ci(1,2); Bias_upper_bound = ci(2,2);
       
        %%% for Radiance
%         plot(L8_band2, D_band1,'.',...
%         [min(L8_band2),max(L8_band2)],ci(1,1)*[min(L8_band2),max(L8_band2)]+ci(1,2),'-b',...
%         [min(L8_band2),max(L8_band2)],ci(2,1)*[min(L8_band2),max(L8_band2)]+ci(2,2),'-g')
%     
%         equation = [num2str(p_slope(1)) '*x + ' num2str(p_slope(2))];
%         tx= strcat('y= ', '  ', equation, ' (fitted red line)');
%         text(20, 240, tx, 'FontSize', 18);
%         tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(p_slope(1)),')',...
%             ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
%         text(20, 230, tx, 'FontSize', 18);
%         tx = strcat('Bias', '(', num2str(p_slope(2)), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
%         text(20, 220, tx, 'FontSize', 18);
%         
%         tx = strcat('Blue line: Lower bound');
%         text(20, 210, tx, 'FontSize', 18);
%         tx = strcat('Green line: Upper bound');
%         text(20, 200, tx, 'FontSize', 18);
        
        %%% for Reflectance
        plot(L8_band2_ref, D_band1_ref,'.',...
        [min(L8_band2_ref),max(L8_band2_ref)],ci(1,1)*[min(L8_band2_ref),max(L8_band2_ref)]+ci(1,2),...
        [min(L8_band2_ref),max(L8_band2_ref)],ci(2,1)*[min(L8_band2_ref),max(L8_band2_ref)]+ci(2,2))
    
        equation = [num2str(p_slope(1)) '*x + ' num2str(p_slope(2))];
        tx= strcat('y= ', '  ', equation, ' (fitted line)');
        text(0.02, 0.55, tx, 'FontSize', 18);
        tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(p_slope(1)),')',...
            ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
        text(0.02, 0.53, tx, 'FontSize', 18);
        tx = strcat('Bias', '(', num2str(p_slope(2)), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
        text(0.02, 0.51, tx, 'FontSize', 18);
        tx = strcat('Blue line: Lower bound');
        text(0.02, 0.49, tx, 'FontSize', 18);
        tx = strcat('Green line: Upper bound');
        text(0.02, 0.47, tx, 'FontSize', 18);

      hold off

     elseif band==2
        % band=2
        L8_band3 = allmean_rad_L8(:,band); D_band2 = allmean_rad_D(:,band);
        L8_band3_ref = allmean_ref_L8(:,band); D_band2_ref = allmean_ref_D(:,band);
        
        % Standard Deviation
        L8_band3_SD = all_SD_rad_L8(:,band); D_band2_SD = all_SD_rad_D(:,band); %Radiance
        L8_band3_ref_SD = all_SD_ref_L8(:,band); D_band2_ref_SD = all_SD_ref_D(:,band); %Reflectance
        
        % Relative Standard Deviation
        L8_band3_RSD = all_RSD_rad_L8(:,band); D_band2_RSD = all_RSD_rad_D(:,band); %Radiance
        L8_band3_ref_RSD = all_RSD_ref_L8(:,band); D_band2_ref_RSD = all_RSD_ref_D(:,band); %Reflectance
        
        sz=50; c=1; cal=4 ; 

        for l=1:15
            
            %%% Radiance
%             p_slope = polyfit(L8_band3(c:cal),D_band2(c:cal),1);
%             p_slope = round(p_slope, 5);
%             fitline_ind = fit(L8_band3(c:cal),D_band2(c:cal), 'poly1');
            
            %%% Reflectance
            fitline_ind = fit(L8_band3_ref(c:cal),D_band2_ref(c:cal), 'poly1');
            p_slope = polyfit(L8_band3_ref(c:cal), D_band2_ref(c:cal),1);
            p_slope = round(p_slope, 5);
            
             % Display evaluated equation y = m*x + b
            equation = [num2str(p_slope(1)) '*x + ' num2str(p_slope(2))];
            slope_band2(l)= p_slope(1); bias_band2(l)=p_slope(2);
            
            G_and_B = coeffvalues(fitline_ind);
            Gain_band2(l) = G_and_B(1,1);
            Bias_band2(l) = G_and_B(1,2);
            ci = confint(fitline_ind);
            ci = round(ci, 5);
            Gain_ci_band2= [ci(1,1) ci(2,1)]; 
            Bias_ci_band2= [ci(1,2) ci(2,2)];
            All_Gain.Gain_ci_band2(l,:) = Gain_ci_band2;
            All_Bias.Bias_ci_band2(l,:) = Bias_ci_band2;
            
             figure(2);
            % Evaluate fit equation using polyval
            %f = polyval(p_slope, L8_band3(c:cal));
%             scatter(L8_band3(c:cal), D_band2(c:cal), sz, marker{l})
%             xlim([0 250]); ylim([0 250]);
%             scatter(L8_band3_ref(c:cal), D_band2_ref(c:cal), sz, marker{l})
%             xlim([0 0.6]); ylim([0 0.6]);
           
            %error = D_band2_RSD(c:cal);
           % errorbar(L8_band3(c:cal), D_band2(c:cal), error, marker{l},'LineWidth',0.25)
          
            
            % Radiance
%             error = D_band2_SD(c:cal);
%             errorbar(L8_band3(c:cal), D_band2(c:cal), error, marker{l},'LineWidth',0.25)
%             xlim([0 250]); ylim([0 250]);
            
% %             Reflectance
            error = D_band2_ref_SD(c:cal);
            errorbar(L8_band3_ref(c:cal), D_band2_ref(c:cal), error, marker{l},'LineWidth',0.25)
            xlim([0 0.6]); ylim([0 0.6]);
            
            hold on
            M = location_str(l);
            c=cal+1;
            cal=cal+4;
%             tx= strcat('y= ', '  ', equation, ' (', char(M), ')');
%             text(10, t, tx, 'FontSize', 12);
%             t= t + 7;
%             fprintf('%s\n', tx)

        end
        
        legend(leg_location,'FontSize',16, 'Location','southeast');
        
        % Radiance
%         [p_slope, S] = polyfit(L8_band3, D_band2,1);
%         p_slope = round(p_slope, 5);
%         f = polyval(p_slope, L8_band3); plot(L8_band3, f, '-r'); hold on
%         fitline = fit(L8_band3, D_band2, 'poly1');
        
        % Reflectance
        [p_slope, S] = polyfit(L8_band3_ref, D_band2_ref,1);
        p_slope = round(p_slope, 5);
        fitline = fit(L8_band3_ref, D_band2_ref, 'poly1');
        f = polyval(p_slope, L8_band3_ref); plot(L8_band3_ref, f); hold on
          
        % Confidence bound
        G_and_B = coeffvalues(fitline);
        Gain_band1 = G_and_B(1,1);
        Bias_band1 = G_and_B(1,2);
        ci = confint(fitline);
        ci = round(ci, 5);
        Gain_lower_bound = ci(1,1); Gain_upper_bound = ci(2,1);
        Bias_lower_bound = ci(1,2); Bias_upper_bound = ci(2,2);
       
        %%% for Radiance
%         plot(L8_band3, D_band2,'.',...
%         [min(L8_band3),max(L8_band3)],ci(1,1)*[min(L8_band3),max(L8_band3)]+ci(1,2), '-b',...
%         [min(L8_band3),max(L8_band3)],ci(2,1)*[min(L8_band3),max(L8_band3)]+ci(2,2), '-g')
%     
%         equation = [num2str(p_slope(1)) '*x + ' num2str(p_slope(2))];
%         tx= strcat('y= ', '  ', equation, ' (fitted red line)');
%         text(20, 240, tx, 'FontSize', 18);
%         tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(p_slope(1)),')',...
%             ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
%         text(20, 230, tx, 'FontSize', 18);
%         tx = strcat('Bias', '(', num2str(p_slope(2)), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
%         text(20, 220, tx, 'FontSize', 18);
%         tx = strcat('Blue line: Lower bound');
%         text(20, 210, tx, 'FontSize', 18);
%         tx = strcat('Green line: Upper bound');
%         text(20, 200, tx, 'FontSize', 18);
%         
        %%% for Reflectance
        plot(L8_band3_ref, D_band2_ref,'.',...
        [min(L8_band3_ref),max(L8_band3_ref)],ci(1,1)*[min(L8_band3_ref),max(L8_band3_ref)]+ci(1,2),...
        [min(L8_band3_ref),max(L8_band3_ref)],ci(2,1)*[min(L8_band3_ref),max(L8_band3_ref)]+ci(2,2))
        equation = [num2str(p_slope(1)) '*x + ' num2str(p_slope(2))];
        tx= strcat('y= ', '  ', equation, ' (fitted line)');
        text(0.02, 0.55, tx, 'FontSize', 18);
        tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(p_slope(1)),')',...
            ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
        text(0.02, 0.53, tx, 'FontSize', 18);
        tx = strcat('Bias', '(', num2str(p_slope(2)), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
        text(0.02, 0.51, tx, 'FontSize', 18);
        tx = strcat('Blue line: Lower bound');
        text(0.02, 0.49, tx, 'FontSize', 18);
        tx = strcat('Green line: Upper bound');
        text(0.02, 0.47, tx, 'FontSize', 18);

        hold off
        
     elseif band==3
        % band = 3
        L8_band4 = allmean_rad_L8(:,band); D_band3 = allmean_rad_D(:,band);
        L8_band4_ref = allmean_ref_L8(:,band); D_band3_ref = allmean_ref_D(:,band);
        
        % Standard Deviation
        L8_band4_SD = all_SD_rad_L8(:,band); D_band3_SD = all_SD_rad_D(:,band); %Radiance
        L8_band4_ref_SD = all_SD_ref_L8(:,band); D_band3_ref_SD = all_SD_ref_D(:,band); %Reflectance
        
        % Relative Standard Deviation
        L8_band4_RSD = all_RSD_rad_L8(:,band); D_band3_RSD = all_RSD_rad_D(:,band); %Radiance
        L8_band4_ref_RSD = all_RSD_ref_L8(:,band); D_band3_ref_RSD = all_RSD_ref_D(:,band); %Reflectance
      
        sz=50; c=1;cal=4 ; 
        
        for l=1:15

            %%% Radiance 
%             p_slope = polyfit(L8_band4(c:cal),D_band3(c:cal),1);
%             p_slope = round(p_slope, 5);
%             fitline_ind = fit(L8_band4(c:cal),D_band3(c:cal), 'poly1');
            
            %%% Reflectance
            p_slope = polyfit(L8_band4_ref(c:cal),D_band3_ref(c:cal),1);
            p_slope = round(p_slope, 5);
            fitline_ind = fit(L8_band4_ref(c:cal),D_band3_ref(c:cal), 'poly1');
            
            slope_band3(l)= p_slope(1); bias_band3(l)= p_slope(2);
            
            % Display evaluated equation y = m*x + b
            equation = [num2str(p_slope(1)) '*x + ' num2str(p_slope(2))];
           
            G_and_B = coeffvalues(fitline_ind);
            Gain_band3(l) = G_and_B(1,1);
            Bias_band3(l) = G_and_B(1,2);
            ci = confint(fitline_ind);
            ci = round(ci, 5);
            Gain_ci_band3= [ci(1,1) ci(2,1)]; 
            Bias_ci_band3= [ci(1,2) ci(2,2)];
            All_Gain.Gain_ci_band3(l,:) = Gain_ci_band3;
            All_Bias.Bias_ci_band3(l,:) = Bias_ci_band3;
            
            % Evaluate fit equation using polyval
           % f = polyval(p_slope, L8_band4(c:cal));
            
            figure(3);
%             scatter(L8_band4(c:cal), D_band3(c:cal), sz, marker{l})
%             xlim([0 250]); ylim([0 250]);
%             scatter(L8_band4_ref(c:cal), D_band3_ref(c:cal), sz, marker{l})
%             xlim([0 0.6]); ylim([0 0.6]);
            
            %%% Radiance
%             error = D_band3_SD(c:cal);
%             errorbar(L8_band4(c:cal), D_band3(c:cal), error, marker{l},'LineWidth',0.25)
%             xlim([0 250]); ylim([0 250]);
            
%             %%% Reflectance
            error = D_band3_ref_SD(c:cal);
            errorbar(L8_band4_ref(c:cal), D_band3_ref(c:cal), error, marker{l},'LineWidth',0.25)
            xlim([0 0.6]); ylim([0 0.6]);
            
            hold on
            M = location_str(l);
            c=cal+1;
            cal=cal+4;
%             tx= strcat('y= ', '  ', equation, ' (', char(M), ')');
%             text(10, t, tx, 'FontSize', 12);
%             t= t + 7;
%             fprintf('%s\n', tx)

        end
        
        legend(leg_location,'FontSize',16, 'Location','southeast');
        hold on

        % Radiance
%         p_slope = polyfit(L8_band4,D_band3,1);
%         p_slope = round(p_slope, 5);
%         f = polyval(p_slope, L8_band4); plot(L8_band4, f, '-r'); hold on
%         fitline = fit(L8_band4, D_band3, 'poly1');
        
        % Reflectance
        [p_slope, S] = polyfit(L8_band4_ref,D_band3_ref,1);
        p_slope = round(p_slope, 5);
        f = polyval(p_slope, L8_band4_ref); plot(L8_band4_ref, f);hold on
        fitline = fit(L8_band4_ref, D_band3_ref, 'poly1');
        
        G_and_B = coeffvalues(fitline);
        Gain_band1 = G_and_B(1,1);
        Bias_band1 = G_and_B(1,2);
        ci = confint(fitline);
        ci = round(ci, 5);
        Gain_lower_bound = ci(1,1); Gain_upper_bound = ci(2,1);
        Bias_lower_bound = ci(1,2); Bias_upper_bound = ci(2,2);
        
        %%% for Radiance
%         plot(L8_band4, D_band3,'.',...
%         [min(L8_band4),max(L8_band4)],ci(1,1)*[min(L8_band4),max(L8_band4)]+ci(1,2),'-b',...
%         [min(L8_band4),max(L8_band4)],ci(2,1)*[min(L8_band4),max(L8_band4)] + ci(2,2),'-g')
%         equation = [num2str(p_slope(1)) '*x + ' num2str(p_slope(2))];
%         tx= strcat('y= ', '  ', equation, ' (fitted red line)');
%         text(20, 240, tx, 'FontSize', 18);
%         tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(p_slope(1)),')',...
%             ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
%         text(20, 230, tx, 'FontSize', 18);
%         tx = strcat('Bias', '(', num2str(p_slope(2)), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
%         text(20, 220, tx, 'FontSize', 18);
%         tx = strcat('Blue line: Lower bound');
%         text(20, 210, tx, 'FontSize', 18);
%         tx = strcat('Green line: Upper bound');
%         text(20, 200, tx, 'FontSize', 18);
%         
        %%% for Reflectance
        plot(L8_band4_ref, D_band3_ref,'.',...
        [min(L8_band4_ref),max(L8_band4_ref)],ci(1,1)*[min(L8_band4_ref),max(L8_band4_ref)]+ci(1,2),...
        [min(L8_band4_ref),max(L8_band4_ref)],ci(2,1)*[min(L8_band4_ref),max(L8_band4_ref)]+ci(2,2))
        
        equation = [num2str(p_slope(1)) '*x + ' num2str(p_slope(2))];
        tx= strcat('y= ', '  ', equation, ' (fitted line)');
        text(0.02, 0.55, tx, 'FontSize', 18);
        tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(p_slope(1)),')',...
            ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
        text(0.02, 0.53, tx, 'FontSize', 18);
        tx = strcat('Bias', '(', num2str(p_slope(2)), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
        text(0.02, 0.51, tx, 'FontSize', 18);
        tx = strcat('Blue line: Lower bound');
        text(0.02, 0.49, tx, 'FontSize', 18);
        tx = strcat('Green line: Upper bound');
        text(0.02, 0.47, tx, 'FontSize', 18);
        hold off
%          
     elseif band==4
        % band = 4 
        L8_band5 = allmean_rad_L8(:,band); D_band4 = allmean_rad_D(:,band);
        L8_band5_ref = allmean_ref_L8(:,band); D_band4_ref = allmean_ref_D(:,band);
        
        % Standard Deviation
        L8_band5_SD = all_SD_rad_L8(:,band); D_band4_SD = all_SD_rad_D(:,band); %Radiance
        L8_band5_ref_SD = all_SD_ref_L8(:,band); D_band4_ref_SD = all_SD_ref_D(:,band); %Reflectance
        
        % Relative Standard Deviation
        L8_band5_RSD = all_RSD_rad_L8(:,band); D_band4_RSD = all_RSD_rad_D(:,band); %Radiance
        L8_band5_ref_RSD = all_RSD_ref_L8(:,band); D_band4_ref_RSD = all_RSD_ref_D(:,band); %Reflectance
        
        sz=50; c=1; cal=4 ; 

        for l=1:15

            %%% Radiance 
%             p_slope = polyfit(L8_band5(c:cal),D_band4(c:cal),1);
%             p_slope = round(p_slope, 5);
%             fitline_ind = fit(L8_band5(c:cal),D_band4(c:cal), 'poly1');
            
            %%% Reflectance
            p_slope = polyfit(L8_band5_ref(c:cal),D_band4_ref(c:cal),1);
            p_slope = round(p_slope, 5);
            fitline_ind = fit(L8_band5_ref(c:cal),D_band4_ref(c:cal), 'poly1');
            
            slope_band4(l)= p_slope(1); bias_band4(l)=p_slope(2);
            G_and_B = coeffvalues(fitline_ind);
            Gain_band4(l) = G_and_B(1,1);
            Bias_band4(l) = G_and_B(1,2);
            ci = confint(fitline_ind);
            ci = round(ci, 5);
            Gain_ci_band4= [ci(1,1) ci(2,1)]; 
            Bias_ci_band4= [ci(1,2) ci(2,2)];
            All_Gain.Gain_ci_band4(l,:) = Gain_ci_band4;
            All_Bias.Bias_ci_band4(l,:) = Bias_ci_band4;
            
            % Display evaluated equation y = m*x + b
            equation = [num2str(p_slope(1)) '*x + ' num2str(p_slope(2))];
            
            % Evaluate fit equation using polyval
            %f = polyval(p_slope, L8_band5(c:cal));
            
            figure(4);
%             scatter(L8_band5(c:cal), D_band4(c:cal), sz, marker{l});
%             xlim([0 250]); ylim([0 250]);
%             scatter(L8_band5_ref(c:cal), D_band4_ref(c:cal), sz, marker{l})
%             xlim([0 0.6]); ylim([0 0.6]);
            
            %%% Radiance
%             error = D_band4_SD(c:cal);
%             errorbar(L8_band5(c:cal), D_band4(c:cal), error, marker{l},'LineWidth',0.25)
%             xlim([0 250]); ylim([0 250]);
            
            %%% Reflectance
            error = D_band4_ref_SD(c:cal);
            errorbar(L8_band5_ref(c:cal), D_band4_ref(c:cal), error, marker{l},'LineWidth',0.25)
            xlim([0 0.6]); ylim([0 0.6]);
%             
            hold on
            M = location_str(l);
            c=cal+1;
            cal=cal+4;
%             tx= strcat('y= ', '  ', equation, ' (', char(M), ')');
%             text(10, t, tx, 'FontSize', 12);
%             t= t + 7;
%             fprintf('%s\n', tx)

        end
         
        legend(leg_location,'FontSize',16, 'Location','southeast');
        hold on
        
        %%% Radiance
%         p_slope = polyfit(L8_band5,D_band4,1);
%         f = polyval(p_slope, L8_band5); plot(L8_band5, f, '-r'); hold on
%         fitline = fit(L8_band5, D_band4, 'poly1');
         
        %%% Reflectance
        [p_slope, S] = polyfit(L8_band5_ref, D_band4_ref, 1); hold on
        p_slope = round(p_slope, 5);
        f = polyval(p_slope, L8_band5_ref); plot(L8_band5_ref, f); hold on
        fitline = fit(L8_band5_ref, D_band4_ref, 'poly1');
        
        G_and_B = coeffvalues(fitline);
        Gain_band1 = G_and_B(1,1);
        Bias_band1 = G_and_B(1,2);
        ci = confint(fitline);
        ci = round(ci, 5);
        Gain_lower_bound = ci(1,1); Gain_upper_bound = ci(2,1);
        Bias_lower_bound = ci(1,2); Bias_upper_bound = ci(2,2);
        
        %%% Radiance
%         plot(L8_band5, D_band4,'.',...
%         [min(L8_band5),max(L8_band5)],ci(1,1)*[min(L8_band5),max(L8_band5)]+ci(1,2),'-b',...
%         [min(L8_band5),max(L8_band5)],ci(2,1)*[min(L8_band5),max(L8_band5)]+ci(2,2),'-g')
%         
%         equation = [num2str(p_slope(1)) '*x + ' num2str(p_slope(2))];
%         tx= strcat('y= ', '  ', equation, ' (fitted red line)');
%         text(20, 240, tx, 'FontSize', 18);
%         tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(p_slope(1)),')',...
%             ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
%         text(20, 230, tx, 'FontSize', 18);
%         tx = strcat('Bias', '(', num2str(p_slope(2)), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
%         text(20, 220, tx, 'FontSize', 18);
%         
%         tx = strcat('Blue line: Lower bound');
%         text(20, 210, tx, 'FontSize', 18);
%         tx = strcat('Green line: Upper bound');
%         text(20, 200, tx, 'FontSize', 18);

        %%% Reflectance
        plot(L8_band5_ref, D_band4_ref,'.',...
        [min(L8_band5_ref),max(L8_band5_ref)],ci(1,1)*[min(L8_band5_ref),max(L8_band5_ref)]+ci(1,2),...
        [min(L8_band5_ref),max(L8_band5_ref)],ci(2,1)*[min(L8_band5_ref),max(L8_band5_ref)]+ci(2,2))
        equation = [num2str(p_slope(1)) '*x + ' num2str(p_slope(2))];
        tx= strcat('y= ', '  ', equation, ' (fitted line)');
        text(0.02, 0.55, tx, 'FontSize', 18);
        tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(p_slope(1)),')',...
            ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
        text(0.02, 0.53, tx, 'FontSize', 18);
        tx = strcat('Bias', '(', num2str(p_slope(2)), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
        text(0.02, 0.51, tx, 'FontSize', 18);
        tx = strcat('Blue line: Lower bound');
        text(0.02, 0.49, tx, 'FontSize', 18);
        tx = strcat('Green line: Upper bound');
        text(0.02, 0.47, tx, 'FontSize', 18);
        hold off
        
   end
   
      %%% Radiancce
%     title(strcat('Mean TOA Radiance Comparison of L8 and Dove', {', '}, strcat(band_name{band}, ' Band')));
%     xlabel('Mean TOA Radiance L8 ROI (W/Sr/m^2/{\mum})')
%     ylabel('Mean TOA Radiance 1047 ROI (W/Sr/m^2/{\mum})')
    
      %%% Reflectance
    title(strcat('Mean TOA Reflectance Comparison of L8 and Dove', {', '}, strcat(band_name{band}, ' Band')));
    xlabel('Mean TOA Reflectance L8 ROI')
    ylabel('Mean TOA Reflectance 1047 ROI')
    
    hold on
    grid on
    grid minor
    ax  = gca;
    ax.FontSize = 24;
    ax.GridColor = 'k';
%   ax.GridAlpha = 0.8;
    ax.MinorGridColor = 'k';
    %ax.FontName = 'Times New Roman';
 end

 
%%  slope vs date
band_name = {'BLUE','GREEN' ,'RED' ,'NIR' };
band_colors={'c','g','r','m'};

 for band=1:4
     
     if band==1
        %figure; 
        subplot(2,2,1);
        
        plot(DoY, slope_band1, '.', 'markers', 28, 'color', band_colors{band});ylim([0 2]);
        %      sorted=sortrows([DoY' slope_band1']);
        %      sorted_x = sorted(:,1);
        %      sorted_y = sorted(:,2);
        %      subplot(2,2,1)
        %     % figure;
        %      plot(sorted_x, sorted_y, '^', 'color', band_colors{band} ); ylim([0 2]);
        %      hold on
        %      line(sorted_x, sorted_y,'LineWidth',1)

        p_slope = polyfit(DoY, slope_band1, 1);
        p_slope = round(p_slope, 5);
        fitline = fit(DoY', slope_band1', 'poly1');
        ci = confint(fitline);
        ci =round(ci, 5);
%         G_and_B = coeffvalues(fitline);
%         Gain_band1(l) = G_and_B(1,1);
%         Bias_band1(l) = G_and_B(1,2);
%         ci = confint(fitline);
%         Gain_ci_band1= [ci(1,1) ci(2,1)]; 
%         Bias_ci_band1= [ci(1,2) ci(2,2)];
%         All_Gain.Gain_ci_band1(l,:) = Gain_ci_band1;
%         All_Bias.Bias_ci_band1(l,:) = Bias_ci_band1;
        
        hold on
        error = All_Gain.Gain_ci_band1(:,2) -  All_Gain.Gain_ci_band1(:,1);
        errorbar(DoY, slope_band1, error/2, 'oc','LineWidth',1); ylim([0 2]);
       % errorbar(DoY, slope_band1, All_Gain.Gain_ci_band1(:,1),All_Gain.Gain_ci_band1(:,2))
        f = polyval(p_slope, DoY);
        hold on
        plot(DoY, f, 'r','LineWidth',1)
        slope = round(p_slope(1), 5);
        bias = round(p_slope(2), 5);
        equation = [num2str(slope) '*x + ' num2str(bias)];
        tx= strcat('y= ', '  ', equation, ' (', 'fitted red line', ')');
        text(245, 1.9, tx, 'FontSize', 10); 
        
        tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(slope),')',...
            ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
        text(245, 1.8, tx, 'FontSize', 10); 
        tx = strcat('Bias', '(', num2str(bias), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
        text(245, 1.7, tx, 'FontSize', 10); 
        
        hold off
        
     elseif band == 2
         %figure; 
         subplot(2,2,2);
         plot(DoY, slope_band2, '.', 'markers', 28,'color', band_colors{band});ylim([0 2]);
         p_slope = polyfit(DoY, slope_band2, 1);
         fitline = fit(DoY', slope_band2', 'poly1');
         ci = confint(fitline);
         ci =round(ci, 5);
         
         hold on
         error = All_Gain.Gain_ci_band2(:,2) -  All_Gain.Gain_ci_band2(:,1);
         errorbar(DoY, slope_band2, error/2, 'og','LineWidth',1)
        
         f = polyval(p_slope, DoY);
         
         slope = round(p_slope(1), 5);
         bias = round(p_slope(2), 5);
         
         hold on
         plot(DoY, f, 'r','LineWidth',1)

         equation = [num2str(slope) '*x + ' num2str(bias)];
         tx= strcat('y= ', '  ', equation, ' (', 'fitted red line', ')');
         text(245, 1.9, tx, 'FontSize', 10); 
        
         tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(slope),')',...
            ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
         text(245, 1.8, tx, 'FontSize', 10); 
         tx = strcat('Bias', '(', num2str(bias), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
         text(245, 1.7, tx, 'FontSize', 10); 
         hold off 
         
     elseif band == 3
         %figure; 
         subplot(2,2,3);
         plot(DoY, slope_band3, '.', 'markers', 28,'color', band_colors{band});ylim([0 2]);
         p_slope = polyfit(DoY, slope_band3, 1);
         
         fitline = fit(DoY', slope_band3', 'poly1');
         ci = confint(fitline);
         ci =round(ci, 5);
         
         slope = round(p_slope(1), 5);
         bias = round(p_slope(2), 5);
         
         hold on
         error = All_Gain.Gain_ci_band3(:,2) -  All_Gain.Gain_ci_band3(:,1);
         errorbar(DoY, slope_band3, error/2, 'or','LineWidth',1)
         
         f = polyval(p_slope, DoY);
         hold on
         plot(DoY, f, 'r','LineWidth',1)
      
         equation = [num2str(slope) '*x + ' num2str(bias)];
         tx= strcat('y= ', '  ', equation, ' (', 'fitted red line', ')');
         text(245, 1.9, tx, 'FontSize', 10); 
        
        tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(slope),')',...
            ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
        text(245, 1.8, tx, 'FontSize', 10); 
        tx = strcat('Bias', '(', num2str(bias), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
        text(245, 1.7, tx, 'FontSize', 10); 
        hold off
         
     elseif band == 4
         %figure; 
         subplot(2,2,4);
         plot(DoY, slope_band4, '.', 'markers', 28, 'color', band_colors{band});ylim([0 2]);
         p_slope = polyfit(DoY, slope_band4, 1);
         
         fitline = fit(DoY', slope_band4', 'poly1');
         ci = confint(fitline);
         ci =round(ci, 5);
         
         slope = round(p_slope(1), 5);
         bias = round(p_slope(2), 5);
         
         hold on
         error = All_Gain.Gain_ci_band4(:,2) -  All_Gain.Gain_ci_band4(:,1);
         errorbar(DoY, slope_band4, error/2, 'om','LineWidth',1)
        
         f = polyval(p_slope, DoY);
         hold on
         plot(DoY, f, 'r','LineWidth',1)

         equation = [num2str(slope) '*x + ' num2str(bias)];
         tx= strcat('y= ', '  ', equation, ' (', 'fitted red line', ')');
         text(245, 1.9, tx, 'FontSize', 10); 
        
         tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(slope),')',...
            ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
         text(245, 1.8, tx, 'FontSize', 10); 
         tx = strcat('Bias', '(', num2str(bias), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
         text(245, 1.7, tx, 'FontSize', 10); 
         hold off
         
     end 
 
    title(strcat('Gain vs. Day of Year', {', '}, strcat(band_name{band}, ' Band')));
    xlabel('Day of Year')
    ylabel('Gain_{(1047/Landsat 8)}')

    hold on
    grid on
    grid minor
    ax  = gca;
    ax.FontSize = 18;
    ax.GridColor = 'k';
%   ax.GridAlpha = 0.8;
    ax.MinorGridColor = 'k';
    %ax.FontName = 'Times New Roman';

 end

%% bias vs date plot
band_name = {'BLUE','GREEN' ,'RED' ,'NIR' };
band_colors={'c','g','r','m'};

for band = 1:4
    
    if band==1
         %figure
         subplot(2,2,1)
         plot(DoY, bias_band1, '.', 'markers', 28,'color', band_colors{band});
         %ylim([-50 50]);
         ylim([-0.1 0.1]);

         p_slope = polyfit(DoY, bias_band1, 1);
         p_slope = round(p_slope, 5);
         
         fitline = fit(DoY', slope_band1', 'poly1');
         ci = confint(fitline);
         ci =round(ci, 5);
         
         hold on
         error = All_Bias.Bias_ci_band1(:,2) -  abs(All_Bias.Bias_ci_band1(:,1));
         errorbar(DoY, bias_band1, error/2, '.c','LineWidth',1);
         
         f = polyval(p_slope, DoY);
         hold on
         plot(DoY, f, 'r','LineWidth',1)

         equation = [num2str(p_slope(1)) '*x + ' num2str(p_slope(2))];
         
%          % For Radiance
%          tx= strcat('y= ', '  ', equation, ' (', 'fitted red line', ')');
%          text(245, 45, tx, 'FontSize', 10); 
%          tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(p_slope(1)),')',...
%             ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
%          text(245, 40, tx, 'FontSize', 10); %for Radiance 
%          tx = strcat('Bias', '(', num2str(p_slope(2)), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
%          text(245, 35, tx, 'FontSize', 10); 

         %For Reflectance
         tx= strcat('y= ', '  ', equation, ' (', 'fitted red line', ')');
         text(245, 0.09, tx, 'FontSize', 10); 
         tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(p_slope(1)),')',...
            ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
         text(245, 0.080, tx, 'FontSize', 10); % for Reflectance
         tx = strcat('Bias', '(', num2str(p_slope(2)), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
         text(245, 0.070, tx, 'FontSize', 10); % for Reflectance
         
       %  hold off
         
    %      hold on
    %      fitline = fit(DoY', bias_band1', 'poly1');
    %      plot(fitline, DoY,'r','LineWidth',1)
        % hold off
         
    elseif band == 2
         %figure
         subplot(2,2,2)
         plot(DoY, bias_band2, '.', 'markers', 28,'color', band_colors{band});
         %ylim([-50 50]);
         ylim([-0.1 0.1]);
         p_slope = polyfit(DoY, bias_band2, 1); 
         p_slope = round(p_slope, 5);
         
         fitline = fit(DoY', slope_band2', 'poly1');
         ci = confint(fitline);
         ci =round(ci, 5);
         
         hold on
         error = All_Bias.Bias_ci_band2(:,2) -  abs(All_Bias.Bias_ci_band2(:,1));
         errorbar(DoY, bias_band2, error/2, 'og','LineWidth',1);
         
         f = polyval(p_slope, DoY);
         hold on
         plot(DoY, f, 'r','LineWidth',1)

         equation = [num2str(p_slope(1)) '*x + ' num2str(p_slope(2))];
         
         %For Radiance
%          tx= strcat('y= ', '  ', equation, ' (', 'fitted red line', ')');
%          text(245, 45, tx, 'FontSize', 10); 
%          tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(p_slope(1)),')',...
%             ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
%          text(245, 40, tx, 'FontSize', 10); 
%          tx = strcat('Bias', '(', num2str(p_slope(2)), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
%          text(245, 35, tx, 'FontSize', 10); 
         
         %For Reflectance
         tx= strcat('y= ', '  ', equation, ' (', 'fitted red line', ')');
         text(245, 0.09, tx, 'FontSize', 10); 
         tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(p_slope(1)),')',...
            ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
         text(245, 0.080, tx, 'FontSize', 10); % for Reflectance
         tx = strcat('Bias', '(', num2str(p_slope(2)), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
         text(245, 0.070, tx, 'FontSize', 10); % for Reflectance
         hold off
         
     elseif band == 3
         %figure
         subplot(2,2,3)
         plot(DoY, bias_band3, '.', 'markers', 28,'color', band_colors{band});
         %ylim([-50 50]);
         ylim([-0.1 0.1]);
         
         p_slope = polyfit(DoY, bias_band3, 1);
         p_slope = round(p_slope, 5);
         
         fitline = fit(DoY', slope_band3', 'poly1');
         ci = confint(fitline);  ci =round(ci, 5);
         hold on
         error = All_Bias.Bias_ci_band3(:,2) -  abs(All_Bias.Bias_ci_band3(:,1));
         errorbar(DoY, bias_band3, error/2, 'or','LineWidth',1)
         
         f = polyval(p_slope, DoY);
         hold on
         plot(DoY, f, 'r','LineWidth',1)

         equation = [num2str(p_slope(1)) '*x + ' num2str(p_slope(2))];
         
         %For Radiance
%          tx= strcat('y= ', '  ', equation, ' (', 'fitted red line', ')');
%          text(245, 45, tx, 'FontSize', 10); 
%          tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(p_slope(1)),')',...
%             ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
%          text(245, 40, tx, 'FontSize', 10); 
%          tx = strcat('Bias', '(', num2str(p_slope(2)), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
%          text(245, 35, tx, 'FontSize', 10); 
         
         %For Reflectance
         tx= strcat('y= ', '  ', equation, ' (', 'fitted red line', ')');
         text(245, 0.09, tx, 'FontSize', 10); 
         tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(p_slope(1)),')',...
            ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
         text(245, 0.080, tx, 'FontSize', 10); % for Reflectance
         tx = strcat('Bias', '(', num2str(p_slope(2)), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
         text(245, 0.070, tx, 'FontSize', 10); % for Reflectance
         
         hold off
         
     elseif band == 4
         %figure;
         subplot(2,2,4)
         plot(DoY, bias_band4, '.', 'markers', 28,'color', band_colors{band});
         %ylim([-10 70]);
         ylim([-0.1 0.2]);
         
         p_slope = polyfit(DoY, bias_band4, 1);
         p_slope = round(p_slope, 5);
         fitline = fit(DoY', slope_band4', 'poly1');
         
         ci = confint(fitline);  ci =round(ci, 5);
         hold on
         error = All_Bias.Bias_ci_band4(:,2) -  abs(All_Bias.Bias_ci_band4(:,1));
         errorbar(DoY, bias_band4, error/2, 'om','LineWidth',1)
         
         f = polyval(p_slope, DoY);
         hold on
         plot(DoY, f, 'r','LineWidth',1)

         equation = [num2str(p_slope(1)) '*x + ' num2str(p_slope(2))];
         %For Radaiance
%          tx= strcat('y= ', '  ', equation, ' (', 'fitted red line', ')');
%          text(245, 62, tx, 'FontSize', 10); 
%          tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(p_slope(1)),')',...
%             ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
%          text(245, 58, tx, 'FontSize', 10); 
%          tx = strcat('Bias', '(', num2str(p_slope(2)), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
%          text(245, 54, tx, 'FontSize', 10); 
         
         %For Reflectance
         tx= strcat('y= ', '  ', equation, ' (', 'fitted red line', ')');
         text(245, 0.19, tx, 'FontSize', 10); 
         tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(p_slope(1)),')',...
            ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
         text(245, 0.175, tx, 'FontSize', 10); % for Reflectance
         tx = strcat('Bias', '(', num2str(p_slope(2)), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
         text(245, 0.16, tx, 'FontSize', 10); % for Reflectance
         hold off
    end 

    title(strcat('Bias vs. Day of Year', {', '}, strcat(band_name{band}, ' Band')));
    xlabel('Day of Year')
    ylabel('Bias')

    hold on
    grid on
    grid minor
    ax  = gca;
    ax.FontSize = 18;
    ax.GridColor = 'k';
   %ax.GridAlpha = 0.8;
    ax.MinorGridColor = 'k';
    ax.FontName = 'Times New Roman';

end

%% Pixel-by-Pixel Plot
%Pixel-by-Pixel Comparison of TOA Radiance of L8 and Dove
band_name = {'CA' ,'BLUE','GREEN' ,'RED' ,'NIR' ,'SWIR1' ,'SWIR2'};
band_colors={'b','c','g','r','m','[0.6667 0 0]','k'};

 for band=2:5 %L8 band
     band=4
     if band==2
        figure; plot(L8_TOArad_band2, D_TOArad_band1, 'color', band_colors{band}, 'LineStyle','None','Marker','.','markers', 18)
     elseif band==3
        figure; plot(L8_TOArad_band3, D_TOArad_band2, 'color', band_colors{band}, 'LineStyle','None','Marker','.','markers', 18)
     elseif band==4
        figure; plot(L8_TOArad_band4, D_TOArad_band3, 'color', band_colors{band}, 'LineStyle','None','Marker','.','markers', 18)
    elseif band==5
        figure; plot(L8_TOArad_band5, D_TOArad_band4, 'color', band_colors{band}, 'LineStyle','None','Marker','.','markers', 18)
     end
     
    title(strcat('Pixel-by-Pixel Comparison of L8 and Dove', {', '}, strcat(band_name{band}, ' Band')));
    xlabel('TOA Radiance L8 ROI (W/Sr/m^2/{\mum})')
    ylabel('TOA Radiance 1047 ROI (W/Sr/m^2/{\mum})')
    hold on
    plot([0 100], [0 100], 'k')
    hold on
    grid on
    grid minor
    ax  = gca;
    ax.FontSize = 30;
    ax.GridColor = 'k';
%   ax.GridAlpha = 0.8;
    ax.MinorGridColor = 'k';
    %ax.FontName = 'Times New Roman';
 end


%%
% BLUE band
plot(L8_ROIrad_b3, D_ROIrand_band2, 'co', 'LineStyle','None','Marker','.','markers', 18)
hold on
plot([0 700], [0 700], 'k')
xlabel('TOA Radiance L8 ROI (W/Sr/m^2/{\mum})')
ylabel('TOA Radiance 1047 ROI (W/Sr/m^2/{\mum})')
title('Pixel-by-Pixel Blue Band');
hold on
grid on
grid minor
ax  = gca;
ax.FontSize = 30;
ax.GridColor = 'k';
%   ax.GridAlpha = 0.8;
ax.MinorGridColor = 'k';



% GREEN band
plot(L8_ROIrad_b3, TOArad_D_b2_final, 'go', 'LineStyle','None','Marker','.','markers', 18)
hold on
plot([0:160], [0:160], 'k')
xlabel('TOA Radiance L8 ROI (W/Sr/m^2/{\mum})')
ylabel('TOA Radiance 1047 ROI (W/Sr/m^2/{\mum})')
title('Pixel-by-Pixel Green Band');
hold on
grid on
grid minor
ax  = gca;
ax.FontSize = 30;
ax.GridColor = 'k';
%   ax.GridAlpha = 0.8;
ax.MinorGridColor = 'k';


%RED band
plot(L8_ROIrad_b4, TOArad_D_b3_final, 'ro', 'LineStyle','None','Marker','.','markers', 18)
hold on
plot([0:150], [0:150], 'k')
xlabel('TOA Radiance L8 ROI (W/Sr/m^2/{\mum})')
ylabel('TOA Radiance 1047 ROI (W/Sr/m^2/{\mum})')
title('Pixel-by-Pixel Red Band');
hold on
grid on
grid minor
ax  = gca;
ax.FontSize = 30;
ax.GridColor = 'k';
%   ax.GridAlpha = 0.8;
ax.MinorGridColor = 'k';

% NIR band
plot(L8_ROIrad_b5, TOArad_D_b4_final, 'mo', 'LineStyle','None','Marker','.','markers', 18)
hold on
plot([0:200], [0:200], 'k')
xlabel('TOA Radiance L8 ROI (W/Sr/m^2/{\mum})')
ylabel('TOA Radiance 1047 ROI (W/Sr/m^2/{\mum})')
title('Pixel-by-Pixel NIR Band');
hold on
grid on
grid minor
ax  = gca;
ax.FontSize = 30;
ax.GridColor = 'k';
%   ax.GridAlpha = 0.8;
ax.MinorGridColor = 'k';

%Plotting all the ToA radiances of Landsat 8
bands = {'01' '02' '03' '04' '05' '06' '07'}; %literally would be b
band_name = {'CA' ,'BLUE','GREEN' ,'RED' ,'NIR' ,'SWIR1' ,'SWIR2'};
band_colors={'b','c','g','r','m','[0.6667 0 0]','k'};

%For plotting 
Load_L8=load('L8_Data.mat');
TOArad_plot_L8.p=[Load_L8.TOARad_L8_all.Band1 Load_L8.TOARad_L8_all.Band2 Load_L8.TOARad_L8_all.Band3.... 
                  Load_L8.TOARad_L8_all.Band4 Load_L8.TOARad_L8_all.Band5 Load_L8.TOARad_L8_all.Band6....
                  Load_L8.TOARad_L8_all.Band7];

 for band=1:7
    %TOArad_plot= strcat('TOArad_band', bands{band});
    plot(Load_L8.DeciYear, TOArad_plot_L8.p(:,band),'LineStyle','None','Marker','*','Color',band_colors{band},'markers', 14)
    title('Temporal Trend of Landsat 8 OLI');
    xlabel('Decimal Year')
    xlim ([2013 2020])
    ylim([0 400])
    ylabel('Mean TOA Radiance (W/Sr/m^2/nm)')
%   legend({'Landsat 8' 'Sentinel 2A' }) 
    legend({'CA' 'Blue' 'Green' 'Red' 'NIR' 'SWIR1' 'SWIR2' }) 
    hold on
    
    grid on
    grid minor
    ax  = gca;
    ax.FontSize = 30;
    ax.GridColor = 'k';
%   ax.GridAlpha = 0.8;
    ax.MinorGridColor = 'k';
    ax.FontName = 'Times New Roman';
 end


%Plotting all the ToA radiances of Sentinel 2A
bands = {'01' '02' '03' '04' '05' '06' '07'}; %literally would be b
band_name = {'CA' ,'BLUE','GREEN' ,'RED' ,'NIR' ,'SWIR1' ,'SWIR2'};
band_colors={'b','c','g','r','m','[0.6667 0 0]','k'};

%For plotting 
Load_S2=load('S2_Data.mat');
TOArad_plot_S2.p=[Load_S2.TOARad_S2_all.Band1 Load_S2.TOARad_S2_all.Band2 Load_S2.TOARad_S2_all.Band3.... 
                  Load_S2.TOARad_S2_all.Band4 Load_S2.TOARad_S2_all.Band5 Load_S2.TOARad_S2_all.Band6....
                  Load_S2.TOARad_S2_all.Band7];

 for band=1:7
    plot(Load_S2.DeciYear, TOArad_plot_S2.p(:,band),'LineStyle','None','Marker','*','Color',band_colors{band},'markers', 14)
    title('Temporal Trend of Sentinel 2A MSI');
    xlabel('Decimal Year')
    xlim ([2013 2020])
    ylim([0 500])
    ylabel('Mean TOA Radiance (W/Sr/m^2/nm)')
%   legend({'Landsat 8' 'Sentinel 2A' }) 
    legend({'CA' 'Blue' 'Green' 'Red' 'NIR' 'SWIR1' 'SWIR2' }) 
    hold on
    
    grid on
    grid minor
    ax  = gca;
    ax.FontSize = 30;
    ax.GridColor = 'k';
%   ax.GridAlpha = 0.8;
    ax.MinorGridColor = 'k';
    ax.FontName = 'Times New Roman';
 end
 
 
%Band by Band Comparison of TOA Radiance of S2A and L8
bands = {'01' '02' '03' '04' '05' '06' '07'}; %literally would be b
band_name = {'CA' ,'BLUE','GREEN' ,'RED' ,'NIR' ,'SWIR1' ,'SWIR2'};
band_colors={'b','c','g','r','m','[0.6667 0 0]','k'};

%For plotting 
Load_L8=load('L8_Data.mat');
TOArad_plot_L8.p=[Load_L8.TOARad_L8_all.Band1 Load_L8.TOARad_L8_all.Band2 Load_L8.TOARad_L8_all.Band3.... 
                  Load_L8.TOARad_L8_all.Band4 Load_L8.TOARad_L8_all.Band5 Load_L8.TOARad_L8_all.Band6....
                  Load_L8.TOARad_L8_all.Band7];
              
Load_S2=load('S2_Data.mat');
TOArad_plot_S2.p=[Load_S2.TOARad_S2_all.Band1 Load_S2.TOARad_S2_all.Band2 Load_S2.TOARad_S2_all.Band3.... 
                  Load_S2.TOARad_S2_all.Band4 Load_S2.TOARad_S2_all.Band5 Load_S2.TOARad_S2_all.Band6....
                  Load_S2.TOARad_S2_all.Band7];

 for band=7
    %TOArad_plot= strcat('TOArad_band', bands{band});
    plot(Load_L8.DeciYear, TOArad_plot_L8.p(:,band),'b', Load_S2.DeciYear, TOArad_plot_S2.p(:,band), 'c'....
       ,'LineStyle','None','Marker','*','markers', 18)
    title(strcat('Temporal Trend of Landsat 8 OLI and Sentinel 2A MSI', {', '}, strcat(band_name{band}, ' Band')));
    xlabel('Decimal Year')
    xlim ([2013 2020])
    ylim([0 500])
    ylabel('Mean TOA Radiance (W/Sr/m^2/nm)')
    legend({'Landsat 8' 'Sentinel 2A' }) 
    hold on
    L8Mean=mean(Load_L8.TOARad_L8_all.Band1);
   % S2Mean=mean(Load_S2.TOARad_S2_all.Band1);
%     str = {'L8 Mean and S2 Mean:' L8Mean S2Mean};
%     text(2018,80, str,'FontSize',20)

    grid on
    grid minor
    ax  = gca;
    ax.FontSize = 30;
    ax.GridColor = 'k';
%   ax.GridAlpha = 0.8;
    ax.MinorGridColor = 'k';
    %ax.FontName = 'Times New Roman';
 end
 
 
 



%inificient way of plotting 
% band1: Coastal/Aerosol
%Day of Year
figure (1)
scatter(DoY, ToA1Mean,'filled','k')
subplot(1,2,1)
plot(DoY, ToA1Mean,'*')
title('Temporal Trend of Coastal/Aerosol Band')
xlabel('Day of Year')
ylabel('TOA Radiance')

subplot(1,2,2)
plot(DoY, ToA1Mean_corr,'*')
title('Corrected Temporal Trend of Coastal/Aerosol Band')
xlabel('Day of Year')
ylabel('TOA Radiance Corrected')
plot(DeciYear, ToA1Mean,'*',x,yn,'ro'); 

% Decimal Year
figure (1)
scatter(DeciYear, ToA1Mean,'filled','k')
subplot(1,2,1)
plot(DeciYear, ToA1Mean,'*')
xlim([2013 2020])
ylim([0 200])
title('Temporal Trend of Coastal/Aerosol Band')
xlabel('Decimal Year')
ylabel('TOA Radiance')

subplot(1,2,2)
plot(DeciYear, ToA1Mean_corr,'*')
xlim([2013 2020])
ylim([0 200])
title('Corrected Temporal Trend of Coastal/Aerosol Band')
xlabel('Decimal Year')
ylabel('TOA Radiance Corrected')


% Sentinel 2 Decimal Year
figure (1)
plot(DeciYear, MeanTOAradS2,'*',DeciYear, MeanTOAradL8,'ro')
xlim([2013 2020])
ylim([0 200])
title('Temporal Trend of Coastal/Aerosol Band')
xlabel('Decimal Year')
ylabel('TOA Radiance (W/Sr/m^2)')

subplot(1,2,2)
plot(DeciYear, ToA1Mean_corr,'*')
xlim([2013 2020])
ylim([0 200])
title('Corrected Temporal Trend of Coastal/Aerosol Band')
xlabel('Decimal Year')
ylabel('TOA Radiance Corrected')


% band2: Blue
figure (2)
scatter(DoY, ToA2Mean,'filled','b')
subplot(1,2,1)
plot(DoY, ToA2Mean,'*')
title('Temporal Trend of Blue Band')
xlabel('Day of Year')
ylabel('ToA Radiance')
subplot(1,2,2)
plot(DoY, ToA1Mean_corr,'*')
title('Corrected Temporal Trend of Blue Band')
xlabel('Day of Year')
ylabel('ToA Radiance Corrected')

% Decimal Year
figure (2)
scatter(DeciYear, ToA2Mean,'filled','k')
subplot(1,2,1)
plot(DeciYear, ToA2Mean,'*')
xlim([2013 2020])
ylim([0 200])
title('Temporal Trend of Blue Band')
xlabel('Decimal Year')
ylabel('ToA Radiance')
subplot(1,2,2)
plot(DeciYear, ToA2Mean_corr,'*')
xlim([2013 2020])
ylim([0 200])
title('Corrected Temporal Trend of Blue Band')
xlabel('Decimal Year')
ylabel('ToA Radiance Corrected')

% band3: Green
figure (3)
scatter(DoY, ToA3Mean,'filled','g')
subplot(1,2,1)
plot(DoY, ToA3Mean,'*')
title('Temporal Trend of Green Band')
xlabel('Day of Year')
ylabel('ToA Radiance')
subplot(1,2,2)
plot(DoY, ToA3Mean_corr,'*')
title('Corrected Temporal Trend of Green Band')
xlabel('Day of Year')
ylabel('ToA Radiance Corrected')

% Decimal Year
figure (3)
scatter(DeciYear, ToA3Mean,'filled','k')
subplot(1,2,1)
plot(DeciYear, ToA3Mean,'*')
xlim([2013 2020])
ylim([0 300])
title('Temporal Trend of Green Band')
xlabel('Decimal Year')
ylabel('ToA Radiance')
subplot(1,2,2)
plot(DeciYear, ToA3Mean_corr,'*')
xlim([2013 2020])
ylim([0 300])
title('Corrected Temporal Trend of Green Band')
xlabel('Decimal Year')
ylabel('ToA Radiance Corrected')

% band4: Red
figure (4)
scatter(DoY, ToA4Mean,'filled','r')
subplot(1,2,1)
plot(DoY, ToA4Mean,'*')
title('Temporal Trend of Red Band')
xlabel('Day of Year')
ylabel('ToA Radiance')
subplot(1,2,2)
plot(DoY, ToA4Mean_corr,'*')
title('Corrected Temporal Trend of Red Band')
xlabel('Day of Year')
ylabel('ToA Radiance Corrected')

% Decimal Year
figure (4)
scatter(DeciYear, ToA4Mean,'filled','r')
subplot(1,2,1)
plot(DeciYear, ToA4Mean,'*')
xlim([2013 2020])
ylim([0 300])
title('Temporal Trend of Red Band')
xlabel('Decimal Year')
ylabel('ToA Radiance')
subplot(1,2,2)
plot(DeciYear, ToA4Mean_corr,'*')
xlim([2013 2020])
ylim([0 300])
title('Corrected Temporal Trend of Red Band')
xlabel('Decimal Year')
ylabel('ToA Radiance Corrected')

% band5: NIR
figure (5)
scatter(DoY, ToA5Mean,'filled','k')
subplot(1,2,1)
plot(DoY, ToA5Mean,'*')
xlim([2013 2020])
ylim([0 300])
title('Temporal Trend of  NIR Band')
xlabel('Day of Year')
ylabel('ToA Radiance')
subplot(1,2,2)
plot(DoY, ToA5Mean_corr,'*')
xlim([2013 2020])
ylim([0 300])
title('Corrected Temporal Trend of NIR Band')
xlabel('Day of Year')
ylabel('ToA Radiance Corrected')

% Decimal Year
figure (5)
scatter(DeciYear, ToA5Mean,'filled','k')
subplot(1,2,1)
plot(DeciYear, ToA5Mean,'*')
xlim([2013 2020])
ylim([0 300])
title('Temporal Trend of NIR Band')
xlabel('Decimal Year')
ylabel('ToA Radiance')
subplot(1,2,2)
plot(DeciYear, ToA5Mean_corr,'*')
xlim([2013 2020])
ylim([0 300])
title('Corrected Temporal Trend of NIR Band')
xlabel('Decimal Year')
ylabel('ToA Radiance Corrected')

% band6: SWIR-1
figure (6)
scatter(DoY, ToA6Mean,'filled','c')
subplot(1,2,1)
plot(DoY, ToA6Mean,'*')
title('Temporal Trend of  SWIR-1 Band')
xlabel('Day of Year')
ylabel('ToA Radiance')
subplot(1,2,2)
plot(DoY, ToA6Mean_corr,'*')
title('Corrected Temporal Trend of SWIR-1 Band')
xlabel('Day of Year')
ylabel('ToA Radiance Corrected')

% Decimal Year
figure (6)
scatter(DeciYear, ToA6Mean,'filled','c')
subplot(1,2,1)
plot(DeciYear, ToA6Mean,'*')
xlim([2013 2020])
ylim([0 100])
title('Temporal Trend of SWIR-1 Band')
xlabel('Decimal Year')
ylabel('ToA Radiance')
subplot(1,2,2)
plot(DeciYear, ToA6Mean_corr,'*')
xlim([2013 2020])
ylim([0 100])
title('Corrected Temporal Trend of SWIR-1 Band')
xlabel('Decimal Year')
ylabel('ToA Radiance Corrected')

% band7: SWIR-2
figure (7)
scatter(DoY, ToA7Mean,'filled','m')
subplot(1,2,1)
plot(DoY, ToA7Mean,'*')
title('Temporal Trend of  SWIR-2 Band')
xlabel('Day of Year')
ylabel('ToA Radiance')
subplot(1,2,2)
plot(DoY, ToA7Mean_corr,'*')
title('Corrected Temporal Trend of SWIR-2 Band')
xlabel('Day of Year')
ylabel('ToA Radiance Corrected')

% Decimal Year
figure (7)
scatter(DeciYear, ToA7Mean,'filled','c')
subplot(1,2,1)
plot(DeciYear, ToA7Mean,'*')
xlim([2013 2020])
ylim([0 100])
title('Temporal Trend of SWIR-2 Band')
xlabel('Decimal Year')
ylabel('ToA Radiance')
subplot(1,2,2)
plot(DeciYear, ToA7Mean_corr,'*')
xlim([2013 2020])
ylim([0 100])
title('Corrected Temporal Trend of SWIR-2 Band')
xlabel('Decimal Year')
ylabel('ToA Radiance Corrected')

% pix2map(R1,1200,1200)