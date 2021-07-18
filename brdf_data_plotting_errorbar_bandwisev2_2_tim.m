%%% data processing on 'with_angle_cluster_13 data' 

%%%update 3:18 sept nahid; Save all the plots in specified folder
%%%Upadte 2: 11 sept Nahid;
%1.path rowwise plot added on the data after BRDF correction or before BRDF correction??
%2.Deleted NaN values beforehand
%3

%%% Update: 1 June,18(Nahid): replaced sigma based filter with hampel filter because
%   staright sigma based filter can damage the potential trend 

%%%%%whats new in version2?
%1. Changed the structure variable name mean to mean_st(mean structure) that was conflicting with the built in mean function 
%2. Option to filter bandwise or whole scene in that date 
%3. Export the whole Cluster13 as a matrix data or structure data in single package
%4. Generate mean staistics with 3 digit precision
%5. Added some options to play with the BRDF models


clc
clear all
close all

%%set constants %%% set values above 999 for exclusion%%%manually set
%%filter_applied flag/counter
data_folder='sentinel 2A result';
filter_applied=0;
spatial_unc_threshold=15;%%uncertainty in percentage
sigma_weight=7; %%%1000 doesnt have any effect
hampel_window_each_side=8;% will filter out the values above [mean+sigmaweight*(sigma)] and below [mean+-igmaweight*(sigma)]
filter_bandwise=false; %true=1;flase=0;
export_uncorrected_data=false;
export_brdf_corrected_data=false;
brdf_model_type='quadratic'; % copy the strings with inverted coma from following 
path_row_wise_plot=false;
distribution_plot=false;
save_plot=false;
save_plot_directory='figure_save';

% % 'linear'	         Model contains an intercept and linear terms for each predictor.
% % 'interactions'	     Model contains an intercept, linear terms, and all products of pairs of distinct predictors (no squared terms).
% % 'purequadratic' 	 Model contains an intercept, linear terms, and squared terms.
% % 'quadratic'     	 Model contains an intercept, linear terms, interactions, and squared terms.
% % 'polyijk'            Model is a polynomial with all terms up to degree i in the first predictor, degree j in the second predictor, etc. Use numerals 0 through 9. For example, 'poly2111' has a constant plus all linear and product terms, and also contains terms with predictor 1 squared.

%%%%%%some userdefined models in wilkinson notation
%Morakots composite:       'Y~x2+y2+x1^2+y1^2+x1*y1'       %Quadratic(SZA,SAA)+Linear (SZA,SAA,VZA,VAA)
         %unc:    %[3.56325095166159;3.42938523210500;1.84992616314360;2.36975524918379;1.36948207573524;1.74050952575631;2.60274263441487]
%conventional 4 angle:     'Y~x1+x2+y1+y2'                 %same as  predefined 'linear' 
         %unc:    %[3.61664469839710;3.47722105366083;1.87701392033493;2.38878789487870;1.38952844727900;1.79311086826008;2.73855772530471]
%*Nahid's modified quad 1:  'Y ~ 1 + x1 + x1:y1 + y1:x2 + x2*y2 + x1^2 +y2^2'
         %unc:    %[3.01056678861194;2.90633039601144;1.69127393944887;2.35400075664358;1.35600287090540;1.62011481421454;2.51487558745989]
%Nahid's modified quad 2:  'Y ~1+y1*x2+x2*y2+x1^2+y2^2' %intercept,x1,y1,x2,y2,x2*y1,x2*y2,x1^2,x2^2 
         %unc:    %[3.01435993656005;2.90510407992934;1.69032477931521;2.35647122771785;1.35707272012595;1.61923706502921;2.56491009521204]
%Nahid's modified quad 3:  'Y ~ 1 + x1*y1 + y1*x2 + x2*y2 + x1^2 + y2^2'  %Intercept,x1,y1,x2,y2,x2*y1,x2*y2,x1:y1,x1^2,x2^2
         %unc:    %[2.98597944149083;2.86566724518224;1.68770175115220;2.35296366436887;1.35598675494953;1.61420284468206;2.49829654966329] 
%*Nahid's modified quad 4:  'Y ~ 1 + x1*y1 + y1*x2 + x2*y2 + x1^2 +y1^2+y2^2' 
         %unc:    %[2.98154953491553;2.86344823179986;1.68136191291148;2.34514816919447;1.35418737924083;1.61333646262334;2.47813220802981]
                                                       
%built in 'quadratic'
         %unc:   %[2.97396588557406;2.84655551265445;1.66309001815401;2.34006754732163;1.35139550735785;1.60211984631465;2.43520158997748]

%%% subscript_1~solar angles ; subscript_2~view angles 
% % x1 = sind(SZA) .* cosd(SAA); 
% % y1 = sind(SZA) .* sind(SAA);
% % x2 = sind(VZA) .* cosd(VAA);
% % y2 = sind(VZA) .* sind(VAA);





files=dir(data_folder);
files([1 2])=[];

%for ss=3:length(files)
mean_vec=[];std_vec=[];acq_date_vec=[];acq_date_str=[];sza=[];saa=[];vza=[];vaa=[];
 %uncomment here and the corresponding end statement if needed path row wise
for f=1:length(files)%ss%[3:6 8:18]%7%[6 7]%[3 5]%%length(files) //6 for 181 40//
            load(fullfile(data_folder,files(f).name));
   
            empty_index=[];
            for i=1:length(stat)
                if isempty(stat(i).mean)
                   empty_index=[i empty_index];
                end
                if isnan(stat(i).mean)
                   empty_index=[i empty_index];
                end
                
            end
            %remove the empty index's due to folder skipping or other issues
            stat(empty_index)=[];
            %%%%%%%%%%%%%%%%%%%%
           
            for b=1:7
            for i=1:length(stat)
                mean_vec2(i,b)=stat(i).mean(b);
                std_vec2(i,b)=stat(i).std(b);
                %acq_date_vec2(i)=str2dec_yr(stat(i).acq_date);
            end
            end
            
            
            for i=1:length(stat)
                acq_date_vec2(i)=str2dec_yr(stat(i).acq_date);
                sza_vec2(i,1)=stat(i).SolZenMean;
                saa_vec2(i,1)=stat(i).SolAzmMean;
                vza_vec2(i,1)=stat(i).SenZenMean;
                vaa_vec2(i,1)=stat(i).SenAzmMean;
            end
            
              
            
            
            for i=1:length(stat)
                acq_date_str2(i,1:10)=strcat(stat(i).acq_date(:,1:4),'-',stat(i).acq_date(:,5:6),'-',stat(i).acq_date(:,7:8));
            end
            
            mean_vec2=mean_vec2';
            std_vec2=std_vec2';
      
            %%%%combine
            mean_vec=[mean_vec mean_vec2];
            std_vec=[std_vec std_vec2];
            acq_date_vec=[acq_date_vec acq_date_vec2];
            acq_date_str=[acq_date_str;acq_date_str2];
            
            sza=[sza;sza_vec2];
            saa=[saa;saa_vec2];
            vza=[vza;vza_vec2];
            vaa=[vaa;vaa_vec2];
            
            clear mean_vec2 acq_date_vec2 std_vec2 acq_date_str2 sza_vec2 saa_vec2 vza_vec2 vaa_vec2
            mean_bands=mean_vec';std_bands=std_vec';unc_bands=((std_vec./mean_vec)*100)';acq_date=acq_date_vec';
            
            % % % % %            %%% path_row wise export chunk
%            C13=[mean_bands std_bands unc_bands sza saa vza vaa acq_date];
%            save(files(f).name,'C13')
%            clear mean_vec acq_date_vec std_vec acq_date_str sza_vec saa_vec vza_vec vaa_vec
           %%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%end %end statement for exporting path_row_wise data
%%%%%%%%%%%%%%%%%%export data section

%%%%%%%%%%%%%%%%%%export data section
if export_uncorrected_data
    %%%matrix
% % % % % % col_1-7 =mean of 7 bands starting from band 1 (mean_bands) ; 
% % % % % % col_8-14=standard deviation  of 7 bands starting from band 1 (std_bands) ;
% % % % % % col_15-21=spatial uncertainty of 7 bands starting from band 1 (unc_bands);
% % % % % % col_22=solar zenith angle (sza);
% % % % % % col_23=solar azimuth angle (saa);
% % % % % % col_24 =view zenith angle (vza);
% % % % % % col_25=view azimuth angle (vaa);
% % % % % % col_26=acquisition date in decimal year (acq_date) ;
    mean_bands=mean_vec';std_bands=std_vec';unc_bands=((std_vec./mean_vec)*100)';acq_date=acq_date_vec';    
    C13=[mean_bands std_bands unc_bands sza saa vza vaa acq_date];
   
    save('August24_C13_S2B_all_path_row_combined','C13')
    %save(strcat('Cluster13_matrix','_update_',date),'Cluster13_matrix')
    %%%structure
    Cluster13_structure.acq_date=acq_date;Cluster13_structure.mean_bands=mean_bands;Cluster13_structure.std_bands=std_bands;
    Cluster13_structure.unc_bands=unc_bands;Cluster13_structure.sza=sza;Cluster13_structure.saa=saa;Cluster13_structure.vza=vza;
    Cluster13_structure.vaa=vaa;
   % save(strcat('Cluster13_structure','_update_',date),'Cluster13_structure')
    %%%
    clear mean_bands std_bands unc_bands acq_date C13 Cluster13_structure
        
end
%%%%%%%%%%%%%%%%%%

  

%%%%find the minimum day and maximum day of the year
acq_date_datenum=[];
for i=1:size(acq_date_str,1)
    acq_date_datenum=[acq_date_datenum datenum(acq_date_str(i,:))];
end
earliest_datenum=min(acq_date_datenum);
most_recent_datenum=max(acq_date_datenum);
span_of_acq_day=most_recent_datenum-earliest_datenum+1;
%%%%%%%%%%%%ending here

%%%%%%filter using the uncertainty threshhold
%1st filter
%%%%based on spatial uncertainty threshhold find the index of data points and remove them%can set different threshhold for
%%%%different bands
if spatial_unc_threshold<999
 filter_applied=filter_applied+1;
end
%%%%%%%%%%%%%%%%%%%%uncertainty vector
unc_vec=100*(std_vec./mean_vec);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for b=1:7
    
        acq.(strcat('band',num2str(b)))=acq_date_vec;
        mean_st.(strcat('band',num2str(b)))=mean_vec(b,:);
        std.(strcat('band',num2str(b)))=std_vec(b,:);
        unc.(strcat('band',num2str(b)))=unc_vec(b,:);
        SZA.(strcat('band',num2str(b)))=sza;
        SAA.(strcat('band',num2str(b)))=saa;
        VZA.(strcat('band',num2str(b)))=vza;
        VAA.(strcat('band',num2str(b)))=vaa;
end       

%%%%%%%%%avg spatia uncertainty embed added 28 august 2018
for bbb=1:7
avg_spat_unc.(strcat('band',num2str(bbb)))= nanmean(unc.(strcat('band',num2str(bbb))));
end
%%%%%%%%%



if filter_bandwise
    
    for b=1:7
        outlier_index.(strcat('band',num2str(b)))=find(unc.(strcat('band',num2str(b)))>spatial_unc_threshold);
        
      
        mean_st.(strcat('band',num2str(b)))(outlier_index.(strcat('band',num2str(b))))=[];
        std.(strcat('band',num2str(b)))(outlier_index.(strcat('band',num2str(b))))=[];
        unc.(strcat('band',num2str(b)))(outlier_index.(strcat('band',num2str(b))))=[];
        acq.(strcat('band',num2str(b)))(outlier_index.(strcat('band',num2str(b))))=[];
        
        SZA.(strcat('band',num2str(b)))(outlier_index.(strcat('band',num2str(b))))=[];
        SAA.(strcat('band',num2str(b)))(outlier_index.(strcat('band',num2str(b))))=[];
        VZA.(strcat('band',num2str(b)))(outlier_index.(strcat('band',num2str(b))))=[];
        VAA.(strcat('band',num2str(b)))(outlier_index.(strcat('band',num2str(b))))=[];
    end
    
end              

%%%%%%move from here
if ~filter_bandwise
     for band=1:7
        outlier_index.(strcat('band',num2str(band)))=find(unc.(strcat('band',num2str(band)))>spatial_unc_threshold);
        
      for b=1:7
        mean_st.(strcat('band',num2str(b)))(outlier_index.(strcat('band',num2str(band))))=[];
        std.(strcat('band',num2str(b)))(outlier_index.(strcat('band',num2str(band))))=[];
        unc.(strcat('band',num2str(b)))(outlier_index.(strcat('band',num2str(band))))=[];
        acq.(strcat('band',num2str(b)))(outlier_index.(strcat('band',num2str(band))))=[];
        
        SZA.(strcat('band',num2str(b)))(outlier_index.(strcat('band',num2str(band))))=[];
        SAA.(strcat('band',num2str(b)))(outlier_index.(strcat('band',num2str(band))))=[];
        VZA.(strcat('band',num2str(b)))(outlier_index.(strcat('band',num2str(band))))=[];
        VAA.(strcat('band',num2str(b)))(outlier_index.(strcat('band',num2str(band))))=[];
      end
    end
    
end

spatial_outliers=outlier_index;%get temporal outlier indices
%%%%%%%%% get temporal freq
for b=1:7
temp_freq(b)=size((acq.(strcat('band',num2str(b)))),2)/span_of_acq_day;    
end

%%%%%%move from here ends

%%%%%%%%%%%%%%calculate temporal stats after 1st spatial uncertainty filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temporal_mean=[];
temporal_std=[];
temporal_uncertainty=[];

for b=1:7
temporal_mean(b,1)=nanmean(mean_st.(strcat('band',num2str(b))));
temporal_std(b,1)=nanstd( mean_st.(strcat('band',num2str(b))));
temporal_uncertainty=100*(temporal_std./temporal_mean);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%2nd filter based on 2 sigma
if sigma_weight<999
   filter_applied=filter_applied+1; 
end   

clear outlier_index
for b=1:7
    
        [yyyy,temp_index_matrix,xmedian,xsigma]=hampel(mean_st.(strcat('band',num2str(b))),hampel_window_each_side,sigma_weight);
        outlier_index.(strcat('band',num2str(b)))=find(temp_index_matrix);
        
       
       
        if filter_bandwise
            mean_st.(strcat('band',num2str(b)))(outlier_index.(strcat('band',num2str(b))))=[];
            std.(strcat('band',num2str(b)))(outlier_index.(strcat('band',num2str(b))))=[];
            unc.(strcat('band',num2str(b)))(outlier_index.(strcat('band',num2str(b))))=[];
            acq.(strcat('band',num2str(b)))(outlier_index.(strcat('band',num2str(b))))=[];
            SZA.(strcat('band',num2str(b)))(outlier_index.(strcat('band',num2str(b))))=[];
            SAA.(strcat('band',num2str(b)))(outlier_index.(strcat('band',num2str(b))))=[];
            VZA.(strcat('band',num2str(b)))(outlier_index.(strcat('band',num2str(b))))=[];
            VAA.(strcat('band',num2str(b)))(outlier_index.(strcat('band',num2str(b))))=[]; 
        end
       
       if ~filter_bandwise
           for c=1:7
                mean_st.(strcat('band',num2str(c)))(outlier_index.(strcat('band',num2str(b))))=[];
                std.(strcat('band',num2str(c)))(outlier_index.(strcat('band',num2str(b))))=[];
                unc.(strcat('band',num2str(c)))(outlier_index.(strcat('band',num2str(b))))=[];
                acq.(strcat('band',num2str(c)))(outlier_index.(strcat('band',num2str(b))))=[];
                SZA.(strcat('band',num2str(c)))(outlier_index.(strcat('band',num2str(b))))=[];
                SAA.(strcat('band',num2str(c)))(outlier_index.(strcat('band',num2str(b))))=[];
                VZA.(strcat('band',num2str(c)))(outlier_index.(strcat('band',num2str(b))))=[];
                VAA.(strcat('band',num2str(c)))(outlier_index.(strcat('band',num2str(b))))=[];     
           end
       end
end
temporal_outliers=outlier_index; %%grab for tracking purpose %22 oct 2018
clear outlier_index


%%%
%grab uncorrected mean
mean_uncorr=mean_st;
%%%
%%%%%%%%% get temporal freq
for b=1:7
temp_freq(b)=size((acq.(strcat('band',num2str(b)))),2)/span_of_acq_day;    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%calculate temporal mean and uncertainty again after filter
temporal_mean=[];
temporal_std=[];
temporal_uncertainty=[];

for b=1:7
temporal_mean(b,1)=nanmean(mean_st.(strcat('band',num2str(b))));
temporal_std(b,1)=nanstd( mean_st.(strcat('band',num2str(b))));
temporal_uncertainty=100*(temporal_std./temporal_mean);
end
%%%%%%%%%%%%%%%%
temporal_mean_before_brdf=temporal_mean;
temporal_unc_before_brdf=temporal_uncertainty;



for b=1:7
  
    if b==1
    m_color='b';
    legend_text='CA';
    end
    
    if b==2
         m_color='c';
         legend_text='Blue';
    end
    if b==3
         m_color='g';
         legend_text='Green';
    end
    if b==4
         m_color='r';
         legend_text='Red';
    end
    if b==5
         m_color='m';
         legend_text='NIR';
    end
    if b==6
         m_color=[0.6667 0 0];
         legend_text='SWIR1';
    end
    if b==7
         m_color='k';
         legend_text='SWIR2';
    end
    figure(b)
%     set(gca,'color','none')
   
    errorbar( acq.(strcat('band',num2str(b))), mean_st.(strcat('band',num2str(b))), std.(strcat('band',num2str(b))),'o','MarkerSize',8, 'MarkerEdgeColor','black',...
        'MarkerFaceColor',m_color,'DisplayName',legend_text,'CapSize',7,'LineWidth',1,'Color',m_color)
    ylim([temporal_mean(b,1)-0.13 temporal_mean(b,1)+0.13])
    xlim([(min(acq_date_vec)-0.1) (max(acq_date_vec)+0.1)]);
    
    %hold on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plot customization
    
                xlabel(['Decimal Year'],'FontSize',15,'FontWeight','bold','Color','k')
                ylabel(['Mean TOA reflectance'],'FontSize',15,'FontWeight','bold','Color','k')
                set(gca,'fontsize',16)
                %set(gca,'color','none')
                legend('show')
                grid on
                if spatial_unc_threshold<999 && sigma_weight<999
                title([ 'Cluster13 mean TOA refl  (L8)       ' 'No. of Filter Applied=' num2str(filter_applied) '        F1:Spat Unc Threshold=' num2str(round(spatial_unc_threshold,4))...
                    '        F2:Sigma Weight='   num2str(round(sigma_weight,4))] );
                end
                
                if filter_applied==0
                title([ 'Cluster13 mean TOA refl  (L8)       ' 'No. of Filter Applied=' num2str(filter_applied) ] );
                end
                
                if spatial_unc_threshold<999 && sigma_weight>999
                title([ 'Cluster13 mean TOA refl  (L8)      ' 'No. of Filter Applied=' num2str(filter_applied) '        F1:Spat Unc Threshold=' num2str(round(spatial_unc_threshold,4))]);
                end
                
                if spatial_unc_threshold>999 && sigma_weight<999
                 title([ 'Cluster13 mean TOA refl  (L8)       ' 'No. of Filter Applied=' num2str(filter_applied) '        F1:Sigma Weight='   num2str(round(sigma_weight,4))] );
                end
                
                dim = [0.132 0.615 .3 .3];
                str = ['Temporal Mean: ' strcat('B',num2str(b),'=') num2str(round(temporal_mean(b,1),3))] ;
                t=annotation('textbox',dim,'String',str,'FitBoxToText','on','EdgeColor','none');
                t.FontSize=18;
                t.FontWeight='bold';

                dim = [0.335 0.615 .3 .3];
                str2 = ['Temporal Unc: ' strcat('B',num2str(b),'=') num2str(round(temporal_uncertainty(b,1),2)) '%'] ;
                u=annotation('textbox',dim,'String',str2,'FitBoxToText','on','EdgeColor','none');
                u.FontSize=18;
                u.FontWeight='bold';
                
%                 dim = [0.53 0.615 .3 .3];
%                str3 = ['Temporal Freq: ' strcat('B',num2str(b),'=') num2str(round((1/temp_freq(b)),2)) ' Day(s)'] ;
%                str3 = ['Temporal Freq: ' strcat('B',num2str(b),'=') num2str(round(temp_freq(b),2)) ' acquisition(s)/day'  ' @every ' num2str(round(1/temp_freq(b),2)) ' Days'  ] ;
%                u=annotation('textbox',dim,'String',str3,'FitBoxToText','on','EdgeColor','none');
%                u.FontSize=18;
%                u.FontWeight='bold';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
   
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plot all band together

for b=1:7
  
    if b==1
    m_color='b';
    legend_text='CA';
    end
    
    if b==2
         m_color='c';
         legend_text='Blue';
    end
    if b==3
         m_color='g';
         legend_text='Green';
    end
    if b==4
         m_color='r';
         legend_text='Red';
    end
    if b==5
         m_color='m';
         legend_text='NIR';
    end
    if b==6
         m_color=[0.6667 0 0];
         legend_text='SWIR1';
    end
    if b==7
         m_color='k';
         legend_text='SWIR2';
    end
    figure(8)
%     set(gca,'color','none')
    grid on
    errorbar( acq.(strcat('band',num2str(b))), mean_st.(strcat('band',num2str(b))), std.(strcat('band',num2str(b))),'o','MarkerSize',8, 'MarkerEdgeColor','black',...
        'MarkerFaceColor',m_color,'DisplayName',legend_text,'CapSize',7,'LineWidth',1,'Color',m_color)
    ylim([0.1 0.8])
    %ylim([temporal_mean(1,1)-0.13 temporal_mean(6,1)+0.13])%take band 6 as highest value and band 1 as lowest value
    xlim([(min(acq_date_vec)-0.1) (max(acq_date_vec)+0.1)]);
    hold on

    
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plot customization
                xlabel(['Decimal Year'],'FontSize',15,'FontWeight','bold','Color','k')
                ylabel(['Mean TOA reflectance'],'FontSize',15,'FontWeight','bold','Color','k')
                set(gca,'fontsize',16)
                %set(gca,'color','none')
                %legend('show')
                
                if  spatial_unc_threshold<999 && sigma_weight<999
                title([ 'Cluster13 mean TOA refl  (L8)       ' 'No. of Filter Applied=' num2str(filter_applied) '        F1:Spat Unc Threshold=' num2str(round(spatial_unc_threshold,4))...
                    '        F2:Sigma Weight='   num2str(round(sigma_weight,4))] );
                end
                
                if  filter_applied==0
                title([ 'Cluster13 mean TOA refl  (L8)       ' 'No. of Filter Applied=' num2str(filter_applied)] );
                end
                
                if  spatial_unc_threshold<999 && sigma_weight>999
                title([ 'Cluster13 mean TOA refl  (L8)       ' 'No. of Filter Applied=' num2str(filter_applied) '        F1:Spat Unc Threshold=' num2str(round(spatial_unc_threshold,4))]);
                end
                
                if  spatial_unc_threshold>999 && sigma_weight<999
                 title([ 'Cluster13 mean TOA refl  (L8)     ' 'No. of Filter Applied=' num2str(filter_applied) '        F1:Sigma Weight='   num2str(round(sigma_weight,4))] );
                end
               
                
                
                dim = [0.132 0.625 .3 .3];
                str = ['Mean: ' 'B1=' num2str(round(temporal_mean(1,1),3)) '  B2=' num2str(round(temporal_mean(2,1),3))...
                    '  B3=' num2str(round(temporal_mean(3,1),3)) '  B4=' num2str(round(temporal_mean(4,1),3))...
                    '  B5= ' num2str(round(temporal_mean(5,1),3)) '  B6= ' num2str(round(temporal_mean(6,1),3))...
                    '  B7= ' num2str(round(temporal_mean(7,1),3))] ;
                t=annotation('textbox',dim,'String',str,'FitBoxToText','on','EdgeColor','none');
                t.FontSize=14;
                %t.FontWeight='bold'

                
                dim = [0.507 0.625 .3 .3];
                str2 = ['Unc: ' 'B1=' num2str(round(temporal_uncertainty(1,1),2)) '%' '  B2=' num2str(round(temporal_uncertainty(2,1),2)) '%'...
                    '  B3=' num2str(round(temporal_uncertainty(3,1),2)) '%' '  B4=' num2str(round(temporal_uncertainty(4,1),2)) '%'...
                    '  B5=' num2str(round(temporal_uncertainty(5,1),2)) '%' ' B6=' num2str(round(temporal_uncertainty(6,1),2))  ...
                      '%' '  B7=' num2str(round(temporal_uncertainty(7,1),2)) '%'] ;
                u=annotation('textbox',dim,'String',str2,'FitBoxToText','on','EdgeColor','none');
                u.FontSize=14;
                
%                 dim = [0.132 0.598 .3 .3];
%                 str3 = ['Temporal Freq(Days): ' 'B1=' num2str(round((1/temp_freq(1)),2))  '  B2=' num2str(round((1/temp_freq(2)),2)) ...
%                     '  B3=' num2str(round((1/temp_freq(3)),2))  '  B4=' num2str(round((1/temp_freq(4)),2)) ...
%                     '  B5=' num2str(round((1/temp_freq(5)),2))  ' B6=' num2str(round((1/temp_freq(6)),2))  ...
%                        '  B7=' num2str(round((1/temp_freq(7)),2)) ] ;
               % u=annotation('textbox',dim,'String',str3,'FitBoxToText','on','EdgeColor','none');
                u.FontSize=14;
                
                grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% plot angles vs decimal date
for b=7
    
    figure(9)
%     set(gca,'color','none')
    plot( acq.(strcat('band',num2str(b))), SZA.(strcat('band',num2str(b)))','o','MarkerSize',10, 'MarkerEdgeColor','black',...
        'MarkerFaceColor','r','DisplayName','Solar Zenith','LineWidth',1,'Color',m_color)
    hold on
    plot( acq.(strcat('band',num2str(b))), SAA.(strcat('band',num2str(b)))','o','MarkerSize',10, 'MarkerEdgeColor','black',...
        'MarkerFaceColor','m','DisplayName','Solar Azimuth','LineWidth',1,'Color',m_color)
    hold on
    plot( acq.(strcat('band',num2str(b))), VZA.(strcat('band',num2str(b)))','o','MarkerSize',10, 'MarkerEdgeColor','black',...
        'MarkerFaceColor','g','DisplayName','View Zenith','LineWidth',1,'Color',m_color)
    hold on
    plot( acq.(strcat('band',num2str(b))), VAA.(strcat('band',num2str(b)))','o','MarkerSize',8, 'MarkerEdgeColor','black',...
        'MarkerFaceColor','c','DisplayName','View Azimuth','LineWidth',1,'Color',m_color)
   
    xlabel(['Decimal Year'],'FontSize',15,'FontWeight','bold','Color','k')
    ylabel(['SZA'],'FontSize',15,'FontWeight','bold','Color','k')
    set(gca,'fontsize',16)
    ylim([-180 180])
    xlim([(min(acq_date_vec)-0.1) (max(acq_date_vec)+0.1)]);
    legend show
    xlabel(['Decimal Year'],'FontSize',15,'FontWeight','bold','Color','k')
    ylabel(['4 angles '],'FontSize',15,'FontWeight','bold','Color','k')
    title('Decimal Year vs (Sun and View angles)')
    grid on
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% brdf correction 
%%%%%%%%%%%%model and correction starts
for i=1:7
    
x1 = sind(SZA.(strcat('band',num2str(i)))) .* cosd(SAA.(strcat('band',num2str(i)))); % SZA.(strcat('band',num2str(b)))
y1 = sind(SZA.(strcat('band',num2str(i)))) .* sind(SAA.(strcat('band',num2str(i))));
x2 = sind(VZA.(strcat('band',num2str(i)))) .* cosd(VAA.(strcat('band',num2str(i))));
y2 = sind(VZA.(strcat('band',num2str(i)))) .* sind(VAA.(strcat('band',num2str(i))));

X=[x1 y1 x2 y2];


Y=mean_st.(strcat('band',num2str(i)))';
tbl=table(x1,y1,x2,y2,Y);
mdl{i}=fitlm (tbl,brdf_model_type) ;
PredictedValues=predict(mdl{i},X);%%be careful with order


ReferenceValue = nanmean(mean_st.(strcat('band',num2str(i)))');%% get the temporal mean as reference value
ActualValues = mean_st.(strcat('band',num2str(i)))';
mean_corr.(strcat('band',num2str(i)))=(ActualValues./PredictedValues) *ReferenceValue;
end

%%%%%%%%%%%%model and correction ends

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%plot the brddf corrected data
mean_st=mean_corr;
%%%%%%%move here

%%%%%%%move end







temporal_mean=[];
temporal_std=[];
temporal_uncertainty=[];

for b=1:7
temporal_mean(b,1)=nanmean(mean_st.(strcat('band',num2str(b))));
temporal_std(b,1)=nanstd( mean_st.(strcat('band',num2str(b))));
temporal_uncertainty=100*(temporal_std./temporal_mean);
end


% % % for b=1:7
% % % temporal_mean(b,1)=nanmean(mean.(strcat('band',num2str(b))));
% % % temporal_std(b,1)=nanstd(mean.(strcat('band',num2str(b))));
% % % temporal_uncertainty=100*(temporal_std./temporal_mean);
% % % end



%%%plot brdf corrected bandwise
for b=1:7
  
    if b==1
    m_color='b';
    legend_text='CA';
    end
    
    if b==2
         m_color='c';
         legend_text='Blue';
    end
    if b==3
         m_color='g';
         legend_text='Green';
    end
    if b==4
         m_color='r';
         legend_text='Red';
    end
    if b==5
         m_color='m';
         legend_text='NIR';
    end
    if b==6
         m_color=[0.6667 0 0];
         legend_text='SWIR1';
    end
    if b==7
         m_color='k';
         legend_text='SWIR2';
    end
    figure(b+9)
%     set(gca,'color','none')
   % plot( acq.(strcat('band',num2str(b))), mean.(strcat('band',num2str(b)))','o','MarkerSize',8, 'MarkerEdgeColor','black',...
    %    'MarkerFaceColor',m_color,'DisplayName',legend_text,'LineWidth',1,'Color',m_color)
errorbar( acq.(strcat('band',num2str(b))), mean_st.(strcat('band',num2str(b)))', std.(strcat('band',num2str(b))),'o','MarkerSize',8, 'MarkerEdgeColor','black',...
'MarkerFaceColor',m_color,'DisplayName',legend_text,'CapSize',7,'LineWidth',1,'Color',m_color)
    
    ylim([temporal_mean(b,1)-0.13 temporal_mean(b,1)+0.13])
    xlim([(min(acq_date_vec)-0.1) (max(acq_date_vec)+0.1)]);
    
    %hold on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plot customization
    
                xlabel(['Decimal Year'],'FontSize',15,'FontWeight','bold','Color','k')
                ylabel(['Mean TOA reflectance'],'FontSize',15,'FontWeight','bold','Color','k')
                set(gca,'fontsize',16)
                %set(gca,'color','none')
                legend('show')
                grid on
                if spatial_unc_threshold<999 && sigma_weight<999
                title([ 'Cluster13 mean TOA refl(L8)(BRDF corrected)       ' 'No. of Filter Applied=' num2str(filter_applied) '        F1:Spat Unc Threshold=' num2str(round(spatial_unc_threshold,4))...
                    '        F2:Sigma Weight='   num2str(round(sigma_weight,4))] );
                end
                
                if filter_applied==0
                title([ 'Cluster13 mean TOA refl(L8)(BRDF corrected)       ' 'No. of Filter Applied=' num2str(filter_applied) ] );
                end
                
                if spatial_unc_threshold<999 && sigma_weight>999
                title([ 'Cluster13 mean TOA refl(L8)(BRDF corrected)       ' 'No. of Filter Applied=' num2str(filter_applied) '        F1:Spat Unc Threshold=' num2str(round(spatial_unc_threshold,4))]);
                end
                
                if spatial_unc_threshold>999 && sigma_weight<999
                 title([ 'Cluster13 mean TOA refl(L8)(BRDF corrected)       ' 'No. of Filter Applied=' num2str(filter_applied) '        F1:Sigma Weight='   num2str(round(sigma_weight,4))] );
                end
                
                dim = [0.132 0.615 .3 .3];
                str = ['Temporal Mean: ' strcat('B',num2str(b),'=') num2str(round(temporal_mean(b,1),3))] ;
                t=annotation('textbox',dim,'String',str,'FitBoxToText','on','EdgeColor','none');
                t.FontSize=18;
                t.FontWeight='bold';

                dim = [0.325 0.615 .3 .3];
                str2 = ['Temporal Unc: ' strcat('B',num2str(b),'=') num2str(round(temporal_uncertainty(b,1),2)) '%'] ;
                u=annotation('textbox',dim,'String',str2,'FitBoxToText','on','EdgeColor','none');
                u.FontSize=18;
                u.FontWeight='bold';
                
%                 dim = [0.50 0.615 .3 .3];
%                 str3 = ['Temporal Freq: ' strcat('B',num2str(b),'=') num2str(round((1/temp_freq(b)),2)) ' Day(s)'] ;
%                 u=annotation('textbox',dim,'String',str3,'FitBoxToText','on','EdgeColor','none');
%                 u.FontSize=18;
%                 u.FontWeight='bold';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
   
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%plot all band together brdf corected
for b=1:7
  
    if b==1
    m_color='b';
    legend_text='CA';
    end
    
    if b==2
         m_color='c';
         legend_text='Blue';
    end
    if b==3
         m_color='g';
         legend_text='Green';
    end
    if b==4
         m_color='r';
         legend_text='Red';
    end
    if b==5
         m_color='m';
         legend_text='NIR';
    end
    if b==6
         m_color=[0.6667 0 0];
         legend_text='SWIR1';
    end
    if b==7
         m_color='k';
         legend_text='SWIR2';
    end
    figure(17)
%     set(gca,'color','none')
    grid on
    %plot( acq.(strcat('band',num2str(b))), mean.(strcat('band',num2str(b))),'o','MarkerSize',8, 'MarkerEdgeColor','black',...
     %   'MarkerFaceColor',m_color,'DisplayName',legend_text,'LineWidth',1,'Color',m_color)
    plot( acq.(strcat('band',num2str(b))), mean_st.(strcat('band',num2str(b)))','o','MarkerSize',8, 'MarkerEdgeColor','black',...
        'MarkerFaceColor',m_color,'DisplayName',legend_text,'Color',m_color) %'LineWidth',1,
    ylim([0.1 0.8])
    xlim([(min(acq_date_vec)-0.1) (max(acq_date_vec)+0.1)]);
    hold on

    
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plot customization
                xlabel(['Decimal Year'],'FontSize',15,'FontWeight','bold','Color','k')
                ylabel(['Mean TOA reflectance'],'FontSize',15,'FontWeight','bold','Color','k')
                set(gca,'fontsize',16)
                %set(gca,'color','none')
                %legend('show')
                
                title([ 'Cluster13 mean TOA refl(L8)(BRDF corrected)   ' 'No. of Filter Applied=' num2str(filter_applied) '     F1:Spat Unc Threshold=' num2str(round(spatial_unc_threshold,4))]);%...
                %    '     F2:Sigma Weight='   num2str(round(sigma_weight,4))] );
                
               
                dim = [0.132 0.625 .3 .3];
                str = ['Mean: ' 'B1=' num2str(round(temporal_mean(1,1),3)) '  B2=' num2str(round(temporal_mean(2,1),3))...
                    '  B3=' num2str(round(temporal_mean(3,1),3)) '  B4=' num2str(round(temporal_mean(4,1),3))...
                    '  B5= ' num2str(round(temporal_mean(5,1),3)) '  B6= ' num2str(round(temporal_mean(6,1),3))...
                    '  B7= ' num2str(round(temporal_mean(7,1),3))] ;
                t=annotation('textbox',dim,'String',str,'FitBoxToText','on','EdgeColor','none');
                t.FontSize=14;
                %t.FontWeight='bold'

                
                dim = [0.507 0.625 .3 .3];
                str2 = ['Unc: ' 'B1=' num2str(round(temporal_uncertainty(1,1),2)) '%' '  B2=' num2str(round(temporal_uncertainty(2,1),2)) '%'...
                    '  B3=' num2str(round(temporal_uncertainty(3,1),2)) '%' '  B4=' num2str(round(temporal_uncertainty(4,1),2)) '%'...
                    '  B5=' num2str(round(temporal_uncertainty(5,1),2)) '%' ' B6=' num2str(round(temporal_uncertainty(6,1),2))  ...
                      '%' '  B7=' num2str(round(temporal_uncertainty(7,1),2)) '%'] ;
                u=annotation('textbox',dim,'String',str2,'FitBoxToText','on','EdgeColor','none');
                u.FontSize=14;
                
                dim = [0.132 0.598 .3 .3];
%                 str3 = ['Temporal Freq(Acq/day): ' 'B1=' num2str(round(temp_freq(1),2))  '  B2=' num2str(round(temp_freq(2),2)) ...
%                     '  B3=' num2str(round(temp_freq(3),2))  '  B4=' num2str(round(temp_freq(4),2)) ...
%                     '  B5=' num2str(round(temp_freq(5),2))  ' B6=' num2str(round(temp_freq(6),2))  ...
%                        '  B7=' num2str(round(temp_freq(7),2)) ] ;
%                     str3 = ['Temporal Freq(Days): ' 'B1=' num2str(round((1/temp_freq(1)),2))  '  B2=' num2str(round((1/temp_freq(2)),2)) ...
%                     '  B3=' num2str(round((1/temp_freq(3)),2))  '  B4=' num2str(round((1/temp_freq(4)),2)) ...
%                     '  B5=' num2str(round((1/temp_freq(5)),2))  ' B6=' num2str(round((1/temp_freq(6)),2))  ...
%                        '  B7=' num2str(round((1/temp_freq(7)),2)) ] ;
%                 u=annotation('textbox',dim,'String',str3,'FitBoxToText','on','EdgeColor','none');
%                 u.FontSize=14;
                grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










if path_row_wise_plot
%%%%%%%%%adopt path row wise plotting
%%%%creates sitewise plot by reading 16 path row information from the
%%%%folder, and plots them band by band
%%%%creates sitewise plot by reading 16 path row information from the
%%%%folder, and plots them band by band

clc
clearvars -except data_folder brdf_model_type temporal_mean temporal_uncertainty;

mkr_mat=['o' '+' '*' '.' 'x' 's' 'd' '^' 'v' '>' '<' 'p' 'h' '<' 'p' 'h'];
color16=[rgb('Red');rgb('Hotpink') ;rgb('Orange'); rgb('Yellow') ;rgb('Brown') ;rgb('Green') ;rgb('Blue'); rgb('Purple') ;...
    rgb('Crimson'); rgb('DeepPink'); rgb('OrangeRed') ;rgb('DarkKhaki'); ...
    rgb('SaddleBrown'); rgb('Lime'); rgb('Teal'); rgb('DeepSkyBlue')];
colors = {'b','c','g','r','m',[0.6667 0 0],'k' };
bands={'CA','Blue','Green','Red','NIR','SWIR1','SWIR2'};


%%set constants %%% set values above 1000 for exclusion%%%manually set
%%filter_applied flag/counter
filter_applied=0; 
spatial_unc_threshold=10000;%%uncertainty in percentage
sigma_weight=10000; % will filter out the values above [mean+sigmaweight*(sigma)] and below [mean+-igmaweight*(sigma)]


files=dir(data_folder);
files(1:2)=[];

mean_vec=[];std_vec=[];acq_date_vec=[];acq_date_str=[];sza=[];saa=[];vza=[];vaa=[];
for f=1:length(files)
            load(fullfile(data_folder,files(f).name));
            empty_index=[];
            for i=1:length(stat)
                if isempty(stat(i).mean)
                   empty_index=[i empty_index];
                end
                if isnan(stat(i).mean)
                   empty_index=[i empty_index];
                end
            end
            %remove the empty index's due to folder skipping or other issues
            stat(empty_index)=[];
            %%%%%%%%%%%%%%%%%%%%
            for b=1:7
            for i=1:length(stat)
                mean_vec2(i,b)=stat(i).mean(b);
                std_vec2(i,b)=stat(i).std(b);
                %acq_date_vec2(i)=str2dec_yr(stat(i).acq_date);
            end
            end
            for i=1:length(stat)
                acq_date_vec2(i)=str2dec_yr(stat(i).acq_date);
                sza_vec2(i,1)=stat(i).SolZenMean;
                saa_vec2(i,1)=stat(i).SolAzmMean;
                vza_vec2(i,1)=stat(i).SenZenMean;
                vaa_vec2(i,1)=stat(i).SenAzmMean;
            end
            for i=1:length(stat)
                acq_date_str2(i,1:10)=strcat(stat(i).acq_date(:,1:4),'-',stat(i).acq_date(:,5:6),'-',stat(i).acq_date(:,7:8));
            end
            mean_vec2=mean_vec2';
            std_vec2=std_vec2';
      %%%%%grab from here
      sites(f).mean=mean_vec2;
      sites(f).std=std_vec2;
      sites(f).acq=acq_date_vec2;
      sites(f).SZA=sza_vec2;
      sites(f).VZA=vza_vec2;
      sites(f).SAA=saa_vec2;
      sites(f).VAA=vaa_vec2;
      
      sites(f).path=files(f).name(end-13:end-11);
      sites(f).row=files(f).name(end-5:end-4);
      %%%%%
            %%%%%%combine // wont need for this program; just for reference
            mean_vec=[mean_vec mean_vec2];
            std_vec=[std_vec std_vec2];
            acq_date_vec=[acq_date_vec acq_date_vec2];
            acq_date_str=[acq_date_str;acq_date_str2];
            sza=[sza;sza_vec2];
            saa=[saa;saa_vec2];
            vza=[vza;vza_vec2];
            vaa=[vaa;vaa_vec2];
            %%%%%%
            clear mean_vec2 acq_date_vec2 std_vec2 acq_date_str2 sza_vec2 saa_vec2 vza_vec2 vaa_vec2
end

%%%sitewise filtering %planned to do later

%%%sitewise brdf correction %inccluded in the loop



%%%%% ok 16 sites seperation done
%%%%% now need to plot them seperately 
for band_num=1:7
    figure(band_num+17)
    hold on
            for s=1:16
               %for i=1:7
    
                    x1 = sind(sites(s).SZA) .* cosd(sites(s).SAA); % SZA.(strcat('band',num2str(b)))
                    y1 = sind(sites(s).SZA) .* sind(sites(s).SAA);
                    x2 = sind(sites(s).VZA) .* cosd(sites(s).VAA);
                    y2 = sind(sites(s).VZA) .* sind(sites(s).VAA);

                    X=[x1 y1 x2 y2];


                    Y=sites(s).mean(band_num,:)';
                    tbl=table(x1,y1,x2,y2,Y);
                    mdl{i}=fitlm (tbl,brdf_model_type) ;
                    PredictedValues=predict(mdl{i},X);%%be careful with order


                    ReferenceValue = nanmean(sites(s).mean(band_num,:)');%% get the temporal mean as reference value
                    ActualValues = sites(s).mean(band_num,:)';
                    sites(s).temp_corr=(ActualValues./PredictedValues) *ReferenceValue;
               % end 
                
               plot(sites(s).acq,sites(s).temp_corr','o','MarkerSize',12,'MarkerFaceColor',color16(s,:),'MarkerEdgeColor',color16(s,:))
               %mean
               legend_text{s}=strcat(sites(s).path,'/',sites(s).row);
            end
            
            legend_text(17:21)={'original mean','upper 1 sigma','lower 1 sigma','upper 3 sigma','lower 3 sigma'};
            set(gca,'fontsize',18)
        
end

%%%%% plot the mean and standard deviation line original
o_mean=temporal_mean';
o_unc=temporal_uncertainty'/100;
o_sd=o_mean.*o_unc;
u_bounds=o_mean+o_sd;
l_bounds=o_mean-o_sd;


%%%max min acq
for s=1:16
max_v(s)=max(sites(s).acq);
min_v(s)=min(sites(s).acq);
end
max_v=max(max_v);
min_v=min(min_v);
%%%%



for band_num=1:7
    figure(band_num+17)
    hold on
    line([min_v max_v],[o_mean(band_num) o_mean(band_num)],'LineWidth',5,'Color',colors{band_num})
    
    line([min_v max_v],[o_mean(band_num)+o_sd(band_num) o_mean(band_num)+o_sd(band_num)],'LineWidth',5,'LineStyle','--','Color',colors{band_num});%astd
    line([min_v max_v],[o_mean(band_num)-o_sd(band_num) o_mean(band_num)-o_sd(band_num)],'LineWidth',5,'LineStyle','--','Color',colors{band_num});%1std
    
    line([min_v max_v],[o_mean(band_num)+(3*o_sd(band_num)) o_mean(band_num)+(3*o_sd(band_num))],'LineWidth',5,'LineStyle',':','Color',colors{band_num});%astd
    line([min_v max_v],[o_mean(band_num)-(3*o_sd(band_num)) o_mean(band_num)-(3*o_sd(band_num))],'LineWidth',5,'LineStyle',':','Color',colors{band_num});%1std
    legend(legend_text)
    %xlim([min_v-0.1 max_v+.5])
    %ylim([(nanmedian(mean_vec(band_num,:))-0.07) (nanmedian(mean_vec(band_num,:))+0.07)])
    xlabel(['Decimal Year'],'FontSize',15,'FontWeight','bold','Color','k')
    ylabel(['Mean TOA reflectance'],'FontSize',15,'FontWeight','bold','Color','k')
    title([ 'Cluster13 (L8) Band' num2str(band_num) '(' bands{band_num} ')' ' ' 'mean TOA reflectance(sitewise and not BRDF Corrected)'] );
    legend(legend_text)
    set(gca,'fontsize',18)
    set(gcf,'color','w')
end

end    


if distribution_plot
%%%%histogram plot
close all
bands={'CA','Blue','Green','Red','NIR','SWIR1','SWIR2'};
for b=1:7
figure(b+25)

subplot(1,2,1)
histogram(mean_corr.(strcat('band',num2str(b))),50);
title(strcat('Distribution of TOA reflectance (L8): ',bands(b)))
subplot(1,2,2)
normplot(mean_corr.(strcat('band',num2str(b))))
end

end

%%%%





if export_brdf_corrected_data
    %initiate
    C13=[];
    
    %gather means
    for b=1:7
    C13(:,b)=mean_corr.(strcat('band',num2str(b)));
    end
    
    %gather std
    for b=1:7
    C13(:,b+7)=std.(strcat('band',num2str(b)));
    end
    
% % % %     %gather uncertainty
% % % %     for b=1:7
% % % %     C13(:,b+7+7)=unc.(strcat('band',num2str(b)));
% % % %     end
% % % %     
% % % %     %gather angles
% % % %     C13(:,b+15)=SZA.(strcat('band',num2str(b)));
% % % %     C13(:,b+16)=SAA.(strcat('band',num2str(b)));
% % % %     C13(:,b+17)=VZA.(strcat('band',num2str(b)));
% % % %     C13(:,b+18)=VAA.(strcat('band',num2str(b)));
% % % %     
% % % %     %gather acquistion date
% % % %     C13(:,b+19)=acq.(strcat('band',num2str(b)));
    C13(:,15)=acq.(strcat('band',num2str(1)));
    
end

   
%%%%%%excel export

% % % % clc
% % % % 
% % % % filename='C13_.xlsx' 
% % % % sheet = 1;
% % % % xlRange = 'A2';
% % % % xlswrite(filename,C13,sheet,xlRange)
% % % % 
% % % % bands={'mean CA','mean Blue','mean Green','mean Red','mean NIR','mean SWIR1','mean SWIR2'...
% % % %     'std CA','std Blue','std Green','std Red','std NIR','std SWIR1','std SWIR2','acq Date'};
% % % % xlRange = 'A1';
% % % % xlswrite(filename,bands,sheet,xlRange)
% % % % 
% % % % %%%%%%%
% % % % 
%%%plot the predicted vs normal
% figure(38);plot(acq.band1,PredictedValues(:,1),'o')
% hold on
% figure(38);plot(acq.band1,ActualValues(:,1),'o')
%%%%
if save_plot
%%%save all figures
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
fig_n=0;
for iFig = 1:length(FigList)
    fig_n=fig_n+1;
  FigHandle = FigList(iFig);
  frame_h=get(FigHandle,'JavaFrame');
  set(frame_h,'Maximized',1);
  FigName   = get(FigHandle, 'Name');
  saveas(FigHandle, fullfile(save_plot_directory, FigName,strcat(num2str(fig_n),'.jpg')));
end
%%%
end