% Ploting for all the bands
band_name = {'BLUE','GREEN' ,'RED' ,'NIR' };
band_colors={'c','g','r','m'};
leg ={'ROI 1', 'ROI 2', 'ROI 3',  'ROI 4', 'ROI 5',  'ROI 6', 'ROI 7'};

all_date = {'2018-12-24','2018-12-24', '2019-01-25', '2019-02-11','2019-02-11', '2019-03-14', '2019-01-22', '2019-01-22', '2019-03-02'};           
location_str = string(leg_location);
%location = {'California', 'Wisconsin', 'Indiana',  'Missouri', 'J-SA'};
marker = {'o', 's', 'd', 'h','p',  '^','>', '<', 'v', 'o', 's', 'd', '+','p','x','*','+', '+', '.'};

%%
for band = 1:4
    if band == 1
    %%% Mean band = 1
    L8_band2 = all_mean_rad_L8(:,band); D_band1 = all_mean_rad_DR(:,band); %Radiance
    L8_band2_ref = all_mean_ref_L8(:,band); D_band1_ref = all_mean_ref_DR(:,band); %Reflectance
        
    %%% Standard Deviation
    %L8_band2_SD = all_SD_rad_L7(:,band); L8_band2_ref_SD = all_SD_ref_L7(:,band); 
    D_band1_SD = all_SD_rad_DR(:,band); %Radiance
    D_band1_ref_SD = all_SD_ref_DR(:,band); %Reflectance
    
    %%% Signal to Noise Ratio
    SNR_band1_D = D_band1./D_band1_SD; 
    
    % Reflectance
    SNR_band1_D_ref = D_band1_ref./D_band1_ref_SD;
   
    sz=50; c=1; cal=4; 
    for l=1:6
        figure(1)
        pos1 = [0.05 0.3 0.42 0.5];
        subplot('Position',pos1)
        
        %%% Radiance
%         p_slope = polyfit(L8_band2(c:cal),D_band1(c:cal),1);
%         p_slope = round(p_slope, 5);
%         fitline_ind = fit(L8_band2(c:cal),D_band1(c:cal), 'poly1');
%         %Error bar
%         error = D_band1_SD(c:cal);
%         errorbar(L8_band2(c:cal), D_band1(c:cal), error, marker{l}, 'MarkerSize',6, 'LineWidth',0.25)
%         xlim([0 160]); ylim([0 160]);
        
        % Reflectance
        p_slope = polyfit(L8_band2_ref(c:cal), D_band1_ref(c:cal), 1);
        p_slope = round(p_slope, 5);
        fitline_ind = fit(L8_band2_ref(c:cal), D_band1_ref(c:cal), 'poly1');
        
        % Error bar
        error = D_band1_ref_SD(c:cal);
        errorbar(L8_band2_ref(c:cal), D_band1_ref(c:cal), error, marker{l},'LineWidth',0.25)
        xlim([0 0.3]); ylim([0 0.3]);

        %Display evaluated equation y = m*x + b
        equation = [num2str(p_slope(1)) '*x + ' num2str(p_slope(2))];
        slope_band1(l)= p_slope(1); bias_band1(l)=p_slope(2);
        
        % Gain, bias, CI
        G_and_B = coeffvalues(fitline_ind);
        Gain_band1(l) = G_and_B(1,1);
        Bias_band1(l) = G_and_B(1,2);
        ci = confint(fitline_ind);
        ci = round(ci, 5);
        Gain_ci_band1= [ci(1,1) ci(2,1)]; 
        Bias_ci_band1= [ci(1,2) ci(2,2)];
        All_Gain.Gain_ci_band1(l,:) = Gain_ci_band1;
        All_Bias.Bias_ci_band1(l,:) = Bias_ci_band1;
      
        %%% Day of Year 
        date = all_date(l);
        [doy,fraction] = date2doy(datenum(char(date)));
        DoY(l)=doy+fraction; %day of year

        c=cal+1; cal=cal+4;  
        hold on
    end

        legend(leg_location,'FontSize',10, 'Location','southeast');
        
        %%% Radiance
%         [p_slope, S] = polyfit(L8_band2, D_band1,1); p_slope = round(p_slope, 5);
%         f = polyval(p_slope, L8_band2); plot(L8_band2, f, '-r'); hold on
%         fitline = fit(L8_band2, D_band1, 'poly1');
        
        %%% Reflectance
        [p_slope, S] = polyfit(L8_band2_ref, D_band1_ref,1);
        p_slope = round(p_slope, 5);
        f = polyval(p_slope, L8_band2_ref); 
        plot(L8_band2_ref, f, '-r'); hold on
        fitline = fit(L8_band2_ref, D_band1_ref, 'poly1');
        
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
%         text(5, 153, tx, 'FontSize', 12);
%         tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(p_slope(1)),')',...
%             ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
%         text(5, 147, tx, 'FontSize', 12);
%         tx = strcat('Bias', '(', num2str(p_slope(2)), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
%         text(5, 141, tx, 'FontSize', 12);
%         
%         tx = strcat('Blue line: Lower bound');
%         text(5, 135, tx, 'FontSize', 12);
%         tx = strcat('Green line: Upper bound');
%         text(5, 129, tx, 'FontSize', 12);
        
        %%% for Reflectance
        plot(L8_band2_ref, D_band1_ref,'.',...
        [min(L8_band2_ref),max(L8_band2_ref)],ci(1,1)*[min(L8_band2_ref),max(L8_band2_ref)]+ci(1,2),'-b',...
        [min(L8_band2_ref),max(L8_band2_ref)],ci(2,1)*[min(L8_band2_ref),max(L8_band2_ref)]+ci(2,2),'-g')
    
        equation = [num2str(p_slope(1)) '*x + ' num2str(p_slope(2))];
        tx= strcat('y= ', '  ', equation, ' (fitted line)');
        text(0.01, 0.29, tx, 'FontSize', 12);
        tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(p_slope(1)),')',...
            ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
        text(0.01, 0.28, tx, 'FontSize', 12);
        tx = strcat('Bias', '(', num2str(p_slope(2)), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
        text(0.01, 0.27, tx, 'FontSize', 12);
        tx = strcat('Blue line: Lower bound');
        text(0.01, 0.26, tx, 'FontSize', 12);
        tx = strcat('Green line: Upper bound');
        text(0.01, 0.25, tx, 'FontSize', 12);
      
    hold off
    elseif band == 2
    %%% Mean band = 2
    L8_band3= all_mean_rad_L8(:,band); D_band2 = all_mean_rad_DR(:,band); %Radiance
    L8_band3_ref = all_mean_ref_L8(:,band); D_band2_ref = all_mean_ref_DR(:,band); %Reflectance
        
    %%% Standard Deviation
    %L8_band2_SD = all_SD_rad_L7(:,band); L8_band2_ref_SD = all_SD_ref_L7(:,band); 
    D_band2_SD = all_SD_rad_DR(:,band); %Radiance
    D_band2_ref_SD = all_SD_ref_DR(:,band); %Reflectance
    
    sz=50; c=1; cal=4; 
    for l=1:6
        %figure(band)
        %pos1 = [0.05 0.3 0.4 0.5];
        pos2 = [0.55 0.3 0.42 0.5];
        subplot('Position',pos2)
        
        %%% Radiance
%         p_slope = polyfit(L8_band3(c:cal),D_band2(c:cal),1);
%         p_slope = round(p_slope, 5);
%         fitline_ind = fit(L8_band3(c:cal),D_band2(c:cal), 'poly1');
%         %%Error bar
%         error = D_band2_SD(c:cal);
%         errorbar(L8_band3(c:cal), D_band2(c:cal), error, marker{l},'LineWidth',0.25)
%         xlim([0 200]); ylim([0 200]);
        
        %%% Reflectance
        p_slope = polyfit(L8_band3_ref(c:cal), D_band2_ref(c:cal), 1);
        p_slope = round(p_slope, 5);
        fitline_ind = fit(L8_band3_ref(c:cal), D_band2_ref(c:cal), 'poly1');
        %%%Error bar
        error = D_band2_ref_SD(c:cal);
        errorbar(L8_band3_ref(c:cal), D_band2_ref(c:cal), error, marker{l},'LineWidth',0.25)
        xlim([0 0.4]); ylim([0 0.4]);

        % Display evaluated equation y = m*x + b
        equation = [num2str(p_slope(1)) '*x + ' num2str(p_slope(2))];
        slope_band2(l)= p_slope(1); bias_band2(l)=p_slope(2);
        
        % Gain, bias, CI
        G_and_B = coeffvalues(fitline_ind);
        Gain_band2(l) = G_and_B(1,1);
        Bias_band2(l) = G_and_B(1,2);
        ci = confint(fitline_ind);
        ci = round(ci, 5);
        Gain_ci_band2= [ci(1,1) ci(2,1)]; 
        Bias_ci_band2= [ci(1,2) ci(2,2)];
        All_Gain.Gain_ci_band2(l,:) = Gain_ci_band2;
        All_Bias.Bias_ci_band2(l,:) = Bias_ci_band2;
      
        c=cal+1; cal=cal+4;  
        hold on
    end
        legend(leg_location,'FontSize',10, 'Location','southeast');
        hold on
        
%         %%% Radiance
%         [p_slope, S] = polyfit(L8_band3, D_band2, 1); p_slope = round(p_slope, 5);
%         f = polyval(p_slope, L8_band3); plot(L8_band3, f, '-r'); hold on
%         fitline = fit(L8_band3, D_band2, 'poly1');
        
        %%% Reflectance
        [p_slope, S] = polyfit(L8_band3_ref, D_band2_ref,1);
        p_slope = round(p_slope, 5);
        fitline = fit(L8_band3_ref, D_band2_ref, 'poly1');
        f = polyval(p_slope, L8_band3_ref); plot(L8_band3_ref, f); hold on
%           
        %%% Confidence bound
        G_and_B = coeffvalues(fitline);
        Gain_band2 = G_and_B(1,1);
        Bias_band2 = G_and_B(1,2);
        ci = confint(fitline);
        ci = round(ci, 5);
        Gain_lower_bound = ci(1,1); Gain_upper_bound = ci(2,1);
        Bias_lower_bound = ci(1,2); Bias_upper_bound = ci(2,2);
%        
%         %%% for Radiance
%         plot(L8_band3, D_band2,'.',...
%         [min(L8_band3),max(L8_band3)],ci(1,1)*[min(L8_band3),max(L8_band3)]+ci(1,2),'-b',...
%         [min(L8_band3),max(L8_band3)],ci(2,1)*[min(L8_band3),max(L8_band3)]+ci(2,2),'-g')
%     
%         equation = [num2str(p_slope(1)) '*x + ' num2str(p_slope(2))];
%         tx= strcat('y= ', '  ', equation, ' (fitted red line)');
%         text(5, 190, tx, 'FontSize', 12);
%         tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(p_slope(1)),')',...
%             ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
%         text(5, 180, tx, 'FontSize', 12);
%         tx = strcat('Bias', '(', num2str(p_slope(2)), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
%         text(5, 170, tx, 'FontSize', 12);
%         
%         tx = strcat('Blue line: Lower bound');
%         text(5, 160, tx, 'FontSize', 12);
%         tx = strcat('Green line: Upper bound');
%         text(5, 150, tx, 'FontSize', 12);
        
        %%% for Reflectance
        plot(L8_band3_ref, D_band2_ref,'.',...
        [min(L8_band3_ref),max(L8_band3_ref)],ci(1,1)*[min(L8_band3_ref),max(L8_band3_ref)]+ci(1,2),'-b',...
        [min(L8_band3_ref),max(L8_band3_ref)],ci(2,1)*[min(L8_band3_ref),max(L8_band3_ref)]+ci(2,2),'-g')
    
        equation = [num2str(p_slope(1)) '*x + ' num2str(p_slope(2))];
        tx= strcat('y= ', '  ', equation, ' (fitted line)');
        text(0.02, 0.39, tx, 'FontSize', 12);
        tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(p_slope(1)),')',...
            ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
        text(0.02, 0.375, tx, 'FontSize', 12);
        tx = strcat('Bias', '(', num2str(p_slope(2)), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
        text(0.02, 0.36, tx, 'FontSize', 12);
        tx = strcat('Blue line: Lower bound');
        text(0.02, 0.345, tx, 'FontSize', 12);
        tx = strcat('Green line: Upper bound');
        text(0.02, 0.33, tx, 'FontSize', 12);
      hold off
    
    elseif band == 3
    %%% Mean band =3
    L8_band4 = all_mean_rad_L8(:,band); D_band3 = all_mean_rad_DR(:,band); %Radiance
    L8_band4_ref = all_mean_ref_L8(:,band); D_band3_ref = all_mean_ref_DR(:,band); %Reflectance
        
    %%% Standard Deviation
    %L8_band3_SD = all_SD_rad_L7(:,band); L8_band2_ref_SD = all_SD_ref_L7(:,band); 
    D_band3_SD = all_SD_rad_DR(:,band); %Radiance
    D_band3_ref_SD = all_SD_ref_DR(:,band); %Reflectance
    
    sz=50; c=1; cal=4; 
    for l=1:6
        figure(2)
        pos1 = [0.05 0.3 0.42 0.5];
        subplot('Position',pos1)

        %%% Radiance
%         p_slope = polyfit(L8_band4(c:cal),D_band3(c:cal),1);
%         p_slope = round(p_slope, 5);
%         fitline_ind = fit(L8_band4(c:cal),D_band3(c:cal), 'poly1');
%         %%Error bar
%         error = D_band3_SD(c:cal);
%         errorbar(L8_band4(c:cal), D_band3(c:cal), error, marker{l},'LineWidth',0.25)
%         xlim([0 200]); ylim([0 200]);
        
        %%% Reflectance
        p_slope = polyfit(L8_band4_ref(c:cal), D_band3_ref(c:cal), 1);
        p_slope = round(p_slope, 5);
        fitline_ind = fit(L8_band4_ref(c:cal), D_band3_ref(c:cal), 'poly1');
        %%%Error bar
        error = D_band3_ref_SD(c:cal);
        errorbar(L8_band4_ref(c:cal), D_band3_ref(c:cal), error, marker{l},'LineWidth',0.25)
        xlim([0 0.5]); ylim([0 0.5]);

        % Display evaluated equation y = m*x + b
        equation = [num2str(p_slope(1)) '*x + ' num2str(p_slope(2))];
        slope_band3(l)= p_slope(1); bias_band3(l)=p_slope(2);
        
        % Gain, bias, CI
        G_and_B = coeffvalues(fitline_ind);
        Gain_band3(l) = G_and_B(1,1);
        Bias_band3(l) = G_and_B(1,2);
        ci = confint(fitline_ind);
        ci = round(ci, 5);
        Gain_ci_band3= [ci(1,1) ci(2,1)]; 
        Bias_ci_band3= [ci(1,2) ci(2,2)];
        All_Gain.Gain_ci_band3(l,:) = Gain_ci_band3;
        All_Bias.Bias_ci_band3(l,:) = Bias_ci_band3;
      
        c=cal+1; cal=cal+4;  
        hold on
    end
        legend(leg_location,'FontSize',10, 'Location','southeast');
        hold on
        
        %%% Radiance
%         [p_slope, S] = polyfit(L8_band4, D_band3, 1); p_slope = round(p_slope, 5);
%         f = polyval(p_slope, L8_band4); plot(L8_band4, f, '-r'); hold on
%         fitline = fit(L8_band4, D_band3, 'poly1');
        
        %%% Reflectance
        [p_slope, S] = polyfit(L8_band4_ref, D_band3_ref,1);
        p_slope = round(p_slope, 5);
        fitline = fit(L8_band4_ref, D_band3_ref, 'poly1');
        f = polyval(p_slope, L8_band4_ref); plot(L8_band4_ref, f); hold on
          
        %%% Confidence bound
        G_and_B = coeffvalues(fitline);
        Gain_band3 = G_and_B(1,1);
        Bias_band3 = G_and_B(1,2);
        ci = confint(fitline);
        ci = round(ci, 5);
        Gain_lower_bound = ci(1,1); Gain_upper_bound = ci(2,1);
        Bias_lower_bound = ci(1,2); Bias_upper_bound = ci(2,2);
       
        %%% for Radiance
%         plot(L8_band4 , D_band3,'.',...
%         [min(L8_band4),max(L8_band4)],ci(1,1)*[min(L8_band4),max(L8_band4)]+ci(1,2),'-b',...
%         [min(L8_band4),max(L8_band4)],ci(2,1)*[min(L8_band4),max(L8_band4)]+ci(2,2),'-g')
%     
%         equation = [num2str(p_slope(1)) '*x + ' num2str(p_slope(2))];
%         tx= strcat('y= ', '  ', equation, ' (fitted red line)');
%         text(5, 190, tx, 'FontSize', 12);
%         tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(p_slope(1)),')',...
%             ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
%         text(5, 180, tx, 'FontSize', 12);
%         tx = strcat('Bias', '(', num2str(p_slope(2)), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
%         text(5, 170, tx, 'FontSize', 12);
%         
%         tx = strcat('Blue line: Lower bound');
%         text(5, 160, tx, 'FontSize', 12);
%         tx = strcat('Green line: Upper bound');
%         text(5, 150, tx, 'FontSize', 12);
        
        %%% for Reflectance
        plot(L8_band4_ref, D_band3_ref,'.',...
        [min(L8_band4_ref),max(L8_band4_ref)],ci(1,1)*[min(L8_band4_ref),max(L8_band4_ref)]+ci(1,2),'-b',...
        [min(L8_band4_ref),max(L8_band4_ref)],ci(2,1)*[min(L8_band4_ref),max(L8_band4_ref)]+ci(2,2),'-g')
    
        equation = [num2str(p_slope(1)) '*x + ' num2str(p_slope(2))];
        tx= strcat('y= ', '  ', equation, ' (fitted line)');
        text(0.02, 0.48, tx, 'FontSize', 12);
        tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(p_slope(1)),')',...
            ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
        text(0.02, 0.465, tx, 'FontSize', 12);
        tx = strcat('Bias', '(', num2str(p_slope(2)), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
        text(0.02, 0.450, tx, 'FontSize', 12);
        tx = strcat('Blue line: Lower bound');
        text(0.02, 0.435, tx, 'FontSize', 12);
        tx = strcat('Green line: Upper bound');
        text(0.02, 0.420, tx, 'FontSize', 12);
      hold off
      
    elseif band == 4
    %%% Mean band = 4
    L8_band5 = all_mean_rad_L8(:,band); D_band4 = all_mean_rad_DR(:,band); %Radiance
    L8_band5_ref = all_mean_ref_L8(:,band); D_band4_ref = all_mean_ref_DR(:,band); %Reflectance
        
    %%% Standard Deviation
    %L8_band2_SD = all_SD_rad_L7(:,band); L8_band2_ref_SD = all_SD_ref_L7(:,band); 
    D_band4_SD = all_SD_rad_DR(:,band); %Radiance
    D_band4_ref_SD = all_SD_ref_DR(:,band); %Reflectance
    
    sz=50; c=1; cal=4; 
    for l=1:6
     
        pos2 = [0.55 0.3 0.42 0.5];
        subplot('Position',pos2)
        %%% Radiance
%         p_slope = polyfit(L8_band5(c:cal),D_band4(c:cal),1);
%         p_slope = round(p_slope, 5);
%         fitline_ind = fit(L8_band5(c:cal),D_band4(c:cal), 'poly1');
%         %Error bar
%         error = D_band4_SD(c:cal);
%         errorbar(L8_band5(c:cal), D_band4(c:cal), error, marker{l},'LineWidth',0.25)
%         xlim([0 200]); ylim([0 200]);
        
        %%% Reflectance
        p_slope = polyfit(L8_band5_ref(c:cal), D_band4_ref(c:cal), 1);
        p_slope = round(p_slope, 5);
        fitline_ind = fit(L8_band5_ref(c:cal), D_band4_ref(c:cal), 'poly1');
        %%%Error bar
        error = D_band4_ref_SD(c:cal);
        errorbar(L8_band5_ref(c:cal), D_band4_ref(c:cal), error, marker{l},'LineWidth',0.25)
        xlim([0 0.6]); ylim([0 0.6]);

        % Display evaluated equation y = m*x + b
        equation = [num2str(p_slope(1)) '*x + ' num2str(p_slope(2))];
        slope_band4(l)= p_slope(1); bias_band4(l)=p_slope(2);
        
        % Gain, bias, CI
        G_and_B = coeffvalues(fitline_ind);
        Gain_band4(l) = G_and_B(1,1);
        Bias_band4(l) = G_and_B(1,2);
        ci = confint(fitline_ind);
        ci = round(ci, 5);
        Gain_ci_band4= [ci(1,1) ci(2,1)]; 
        Bias_ci_band4= [ci(1,2) ci(2,2)];
        All_Gain.Gain_ci_band3(l,:) = Gain_ci_band4;
        All_Bias.Bias_ci_band3(l,:) = Bias_ci_band4;
      
        c=cal+1; cal=cal+4;  
        hold on
    end
    
        legend(leg_location,'FontSize',10, 'Location','southeast');
        hold on
        
        %%% Radiance
%         [p_slope, S] = polyfit(L8_band5, D_band4, 1); p_slope = round(p_slope, 5);
%         f = polyval(p_slope, L8_band5); plot(L8_band5, f, '-r'); hold on
%         fitline = fit(L8_band5, D_band4, 'poly1');
        
        %%% Reflectance
        [p_slope, S] = polyfit(L8_band5_ref, D_band4_ref,1);
        p_slope = round(p_slope, 5);
        fitline = fit(L8_band5_ref, D_band4_ref, 'poly1');
        f = polyval(p_slope, L8_band5_ref); plot(L8_band5_ref, f); hold on
%           
        %%% Confidence bound
        G_and_B = coeffvalues(fitline);
        Gain_band4 = G_and_B(1,1);
        Bias_band4 = G_and_B(1,2);
        ci = confint(fitline);
        ci = round(ci, 5);
        Gain_lower_bound = ci(1,1); Gain_upper_bound = ci(2,1);
        Bias_lower_bound = ci(1,2); Bias_upper_bound = ci(2,2);
       
        %%% for Radiance
%         plot(L8_band5, D_band4,'.',...
%         [min(L8_band5),max(L8_band5)],ci(1,1)*[min(L8_band5),max(L8_band5)]+ci(1,2),'-b',...
%         [min(L8_band5),max(L8_band5)],ci(2,1)*[min(L8_band5),max(L8_band5)]+ci(2,2),'-g')
%     
%         equation = [num2str(p_slope(1)) '*x + ' num2str(p_slope(2))];
%         tx= strcat('y= ', '  ', equation, ' (fitted red line)');
%         text(5, 190, tx, 'FontSize', 12);
%         tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(p_slope(1)),')',...
%             ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
%         text(5, 180, tx, 'FontSize', 12);
%         tx = strcat('Bias', '(', num2str(p_slope(2)), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
%         text(5, 170, tx, 'FontSize', 12);
%         
%         tx = strcat('Blue line: Lower bound');
%         text(5, 160, tx, 'FontSize', 12);
%         tx = strcat('Green line: Upper bound');
%         text(5, 150, tx, 'FontSize', 12);
        
        %%% for Reflectance
        plot(L8_band5_ref, D_band4_ref,'.',...
        [min(L8_band5_ref),max(L8_band5_ref)],ci(1,1)*[min(L8_band5_ref),max(L8_band5_ref)]+ci(1,2),'-b',...
        [min(L8_band5_ref),max(L8_band5_ref)],ci(2,1)*[min(L8_band5_ref),max(L8_band5_ref)]+ci(2,2),'-g')
    
        equation = [num2str(p_slope(1)) '*x + ' num2str(p_slope(2))];
        tx= strcat('y= ', '  ', equation, ' (fitted line)');
        text(0.02, 0.58, tx, 'FontSize', 12);
        tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(p_slope(1)),')',...
            ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
        text(0.02, 0.56, tx, 'FontSize', 12);
        tx = strcat('Bias', '(', num2str(p_slope(2)), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
        text(0.02, 0.54, tx, 'FontSize', 12);
        tx = strcat('Blue line: Lower bound');
        text(0.02, 0.52, tx, 'FontSize', 12);
        tx = strcat('Green line: Upper bound');
        text(0.02, 0.50, tx, 'FontSize', 12);
      hold off
      
    end
      
      %%% Radiancce
%     title(strcat('Mean TOA Radiance Comparison of L8 and Dove-R(1058)', {', '}, strcat(band_name{band}, ' Band')));
%     xlabel('Mean TOA Radiance L8 ROI (W/Sr/m^2/{\mum})', 'fontweight', 'bold','fontsize',38)
%     ylabel('Mean TOA Radiance D-R(1058) ROI (W/Sr/m^2/{\mum})', 'fontweight', 'bold','fontsize',28)
    
      %%% Reflectance
    title(strcat('Mean TOA Reflectance Comparison of L7 and Dove-R(1068)', {', '}, strcat(band_name{band}, ' Band')));
    xlabel('Mean TOA Reflectance L7 ROI', 'fontweight', 'bold','fontsize',38)
    ylabel('Mean TOA Reflectance Dove-R(1068) ROI', 'fontweight', 'bold','fontsize',38)
    
    hold on
    plot([0 200], [0 200], 'color', 'k', 'LineStyle', '--','LineWidth', 0.5)
    grid on
    grid minor
    ax  = gca;
    ax.FontSize = 16;
    ax.GridColor = 'k';
%   ax.GridAlpha = 0.8;
    ax.MinorGridColor = 'k';
    %ax.FontName = 'Times New Roman';
end

%% SNR Dove R
 for band=1:4
   if band==1    
        %% band =1
        %%% Mean
        %L8_band2 = all_mean_rad_L8(:,band); 
        %L8_band2_ref = all_mean_ref_L8(:,band); 
        
        D_band1 = all_mean_rad_DR(:,band); %Radiance
        D_band1_ref = all_mean_ref_DR(:,band); %Reflectance
        
        %%% Standard Deviation
        D_band1_SD = all_SD_rad_DR(:,band); %Radiance
        D_band1_ref_SD = all_SD_ref_DR(:,band); %Reflectance
        
        %%% Signal to Noise Ratio
        SNR_band1_D = D_band1./D_band1_SD; 
        %SNR_band1_L7 = L8_band2./L8_band2_SD;
        
        % Reflectance
        SNR_band1_D_ref = D_band1_ref./D_band1_ref_SD;
        %SNR_band1_L7_ref = L8_band2_ref./L8_band2_ref_SD; 
        figure
        sz=50; c=1; cal=4; 
        for l=1:6
            %figure(1)
            subplot(2,2,1)
            % Radiance
            scatter(D_band1(c:cal), SNR_band1_D(c:cal), sz, marker{l},'LineWidth',1)
            xlim([0 250]); ylim([0 150]);

            % Reflectance
            %scatter(D_band1_ref(c:cal), SNR_band1_D_ref(c:cal), sz, marker{l},'LineWidth',1)
            %scatter(L8_band2_ref(c:cal), SNR_band2_L7_ref(c:cal), sz, marker{l},'LineWidth',1)
            %xlim([0 0.5]); ylim([0 120]);
            c=cal+1; cal=cal+4;  
            hold on
        end 
        
        legend(leg_location,'FontSize',10, 'Location','southeast');
        
%         fitline = fit(D_band1, SNR_band1_D, 'poly2');
%         figure, plot(fitline, D_band1, '.')
%         
%         [p_slope, S] = polyfit(D_band1, SNR_band1_D, 1);
%         p_slope = round(p_slope, 5);
%         f = polyval(p_slope, D_band1); plot(D_band1, f); hold on

         D_MeanTOA_SNR = [D_band1  SNR_band1_D];
%         D_MeanTOA_SNR_sorted = sortrows(D_MeanTOA_SNR);
%         ii = 5;
%         for i=1:6:length(D_MeanTOA_SNR_sorted)
%                Envelope_data(i) = max(D_MeanTOA_SNR_sorted(1:ii, 2));
%                ii = ii + 4;
%         end
%         
%         Sub = D_MeanTOA_SNR(D_MeanTOA_SNR(:, 1) > 50 & D_MeanTOA_SNR(:, 1) < 70, :);
%         max(Sub(:,2))
%      [up,dw]=envelope(D_band1, 100,'peak'); 
%      figure, plot(D_band1,[up;dw])
%      figure, plot(D_band1, [SNR_band1_D; up; dw])
%      figure, plot(D_band1, SNR_band1_D, '*'); xlim([0 200]); ylim([0 150]); hold on;
%      
%      plot(D_band1, dw, '.')
%      plot(D_band1, up, '.')
        
    elseif band ==2
        % band=2
        L7_band2 = all_mean_rad_L7(:,band); D_band2 = all_mean_rad_DR(:,band);
        L7_band2_ref = all_mean_ref_L7(:,band); D_band2_ref = all_mean_ref_DR(:,band);
        
        % Standard Deviation
        D_band2_SD = all_SD_rad_DR(:,band); %Radiance
        D_band2_ref_SD = all_SD_ref_DR(:,band); %Reflectance
        
        % Reflectance
        SNR_band2_D = D_band2./D_band2_SD;
        SNR_band2_D_ref = D_band2_ref./D_band2_ref_SD;
        
        sz=50; c=1; cal=4; 
        for l=1:6
            %figure(2)
            subplot(2,2,2)
            scatter(D_band2(c:cal), SNR_band2_D(c:cal), sz, marker{l},'LineWidth',1)
            xlim([0 250]); ylim([0 150]);
            
            % Reflectance
           % scatter(D_band2_ref(c:cal), SNR_band2_D_ref(c:cal), sz, marker{l},'LineWidth',1)
%             scatter(L7_band2_ref(c:cal), SNR_band3_L7_ref(c:cal), sz, marker{l},'LineWidth',1)
%             xlim([0 0.6]); ylim([0 100]);
%             
            c=cal+1; cal=cal+4;  
            hold on
        end 
        %legend(leg_location,'FontSize',10, 'Location','southeast');
        
     elseif band == 3
         % band=3
        L7_band3 = all_mean_rad_L7(:,band); D_band3 = all_mean_rad_DR(:,band);
        L7_band3_ref = all_mean_ref_L7(:,band); D_band3_ref = all_mean_ref_DR(:,band);
        
        % Standard Deviation
        D_band3_SD = all_SD_rad_DR(:,band); %Radiance
        D_band3_ref_SD = all_SD_ref_DR(:,band); %Reflectance

        %%% Signal to Noise Ratio
        SNR_band3_D = D_band3./D_band3_SD; 
        
        % Reflectance
        SNR_band3_D_ref = D_band3_ref./D_band3_ref_SD;
        
        sz=50; c=1; cal=4; 
        for l=1:6
            %figure(3)
            subplot(2,2,3)
            % Radiance
            scatter(D_band3(c:cal), SNR_band3_D(c:cal) , sz, marker{l},'LineWidth',1)
            xlim([0 250]); ylim([0 100]);
            
            % Reflectance
            %scatter(D_band3_ref(c:cal), SNR_band3_D_ref(c:cal), sz, marker{l},'LineWidth',1)
%             scatter(L7_band3_ref(c:cal), SNR_band4_L7_ref(c:cal), sz, marker{l},'LineWidth',1)  
%             xlim([0 0.6]); ylim([0 100]);
            
            c=cal+1; cal=cal+4;  
            hold on
        end 
        %legend(leg_location,'FontSize',10, 'Location','southeast');
        
    elseif band == 4
        % band=4
        L7_band4 = all_mean_rad_L7(:,band); D_band4 = all_mean_rad_DR(:,band);
        L7_band4_ref = all_mean_ref_L7(:,band); D_band4_ref = all_mean_ref_DR(:,band);
        
        % Standard Deviation
        D_band4_SD = all_SD_rad_DR(:,band); %Radiance
        D_band4_ref_SD = all_SD_ref_DR(:,band); %Reflectance
        
        %%% Signal to Noise Ratio
        SNR_band4_D = D_band4./D_band4_SD; 
        
        % Reflectance
        SNR_band4_D_ref = D_band4_ref./D_band4_ref_SD; 
        
        sz=50; c=1; cal=4; 
        for l=1:6
            %figure(4)
            subplot(2,2,4)
            % Radiance
            scatter(D_band4(c:cal), SNR_band4_D(c:cal) , sz, marker{l},'LineWidth',1)
            xlim([0 250]); ylim([0 100]);
            
            % Reflectance
            %scatter(D_band4_ref(c:cal), SNR_band4_D_ref(c:cal), sz, marker{l},'LineWidth',1)
%             scatter(L7_band4_ref(c:cal), SNR_band5_L7_ref(c:cal), sz, marker{l},'LineWidth',1)
%             xlim([0 1]); ylim([0 120]);
            
            c=cal+1; cal=cal+4;  
            hold on
        end
        %legend(leg_location,'FontSize',10, 'Location','southeast');
        
   end
   
   
    %%% Dove   
    % Radiance
    title(strcat('SNR vs. Mean TOA Radiance', {', '}, strcat(band_name{band}, ' Band')));
    xlabel('Mean TOA Radiance Dove-R(1058) ROI (W/Sr/m^2/{\mum})')
    ylabel('Signal-to-Noise Ratio')
    
    % Reflectance
%     title(strcat('SNR vs. Mean TOA Reflectance', {', '}, strcat(band_name{band}, ' Band')));
%     xlabel('Mean TOA Reflectance of 1047 ROI')
%     ylabel('Signal-to-Noise Ratio')
%     
    %%% Landsat 8
    % Radiance
%     title(strcat('SNR vs. Mean TOA Radiance', {', '}, strcat(band_name{band}, ' Band')));
%     xlabel('Mean TOA Radiance L7 ROI (W/Sr/m^2/{\mum})')
%     ylabel('Signal-to-Noise Ratio')
    
    % Reflectance
%     title(strcat('SNR vs. Mean TOA Reflectance', {', '}, strcat(band_name{band}, ' Band')));
%     xlabel('Mean TOA Reflectance of L7 ROI')
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
        text(20, 1.9, tx, 'FontSize', 10); 
        
        tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(slope),')',...
            ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
        text(20, 1.8, tx, 'FontSize', 10); 
        tx = strcat('Bias', '(', num2str(bias), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
        text(20, 1.7, tx, 'FontSize', 10); 
        
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
         text(20, 1.9, tx, 'FontSize', 10); 
        
         tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(slope),')',...
            ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
         text(20, 1.8, tx, 'FontSize', 10); 
         tx = strcat('Bias', '(', num2str(bias), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
         text(20, 1.7, tx, 'FontSize', 10); 
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
         text(20, 1.9, tx, 'FontSize', 10); 
        
        tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(slope),')',...
            ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
        text(20, 1.8, tx, 'FontSize', 10); 
        tx = strcat('Bias', '(', num2str(bias), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
        text(20, 1.7, tx, 'FontSize', 10); 
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
         text(20, 1.9, tx, 'FontSize', 10); 
        
         tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(slope),')',...
            ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
         text(20, 1.8, tx, 'FontSize', 10); 
         tx = strcat('Bias', '(', num2str(bias), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
         text(20, 1.7, tx, 'FontSize', 10); 
         hold off
         
     end 
 
    title(strcat('Gain vs. Day of Year', {', '}, strcat(band_name{band}, ' Band')));
    xlabel('Day of Year')
    ylabel('Gain_{(Dove-R/Landsat 8)}')

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
        % ylim([-0.1 0.1]);

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
         tx= strcat('y= ', '  ', equation, ' (', 'fitted red line', ')');
         text(20, 45, tx, 'FontSize', 10); 
         tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(p_slope(1)),')',...
            ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
         text(20, 40, tx, 'FontSize', 10); %for Radiance 
         tx = strcat('Bias', '(', num2str(p_slope(2)), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
         text(20, 35, tx, 'FontSize', 10); 

         %For Reflectance
%          tx= strcat('y= ', '  ', equation, ' (', 'fitted red line', ')');
%          text(245, 0.09, tx, 'FontSize', 10); 
%          tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(p_slope(1)),')',...
%             ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
%          text(245, 0.080, tx, 'FontSize', 10); % for Reflectance
%          tx = strcat('Bias', '(', num2str(p_slope(2)), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
%          text(245, 0.070, tx, 'FontSize', 10); % for Reflectance
         
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
         %ylim([-0.1 0.1]);
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
         tx= strcat('y= ', '  ', equation, ' (', 'fitted red line', ')');
         text(20, 45, tx, 'FontSize', 10); 
         tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(p_slope(1)),')',...
            ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
         text(20, 40, tx, 'FontSize', 10); 
         tx = strcat('Bias', '(', num2str(p_slope(2)), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
         text(20, 35, tx, 'FontSize', 10); 
         
         %For Reflectance
%          tx= strcat('y= ', '  ', equation, ' (', 'fitted red line', ')');
%          text(20, 0.09, tx, 'FontSize', 10); 
%          tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(p_slope(1)),')',...
%             ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
%          text(20, 0.080, tx, 'FontSize', 10); % for Reflectance
%          tx = strcat('Bias', '(', num2str(p_slope(2)), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
%          text(20, 0.070, tx, 'FontSize', 10); % for Reflectance
%          hold off
         
     elseif band == 3
         %figure
         subplot(2,2,3)
         plot(DoY, bias_band3, '.', 'markers', 28,'color', band_colors{band});
         %ylim([-50 50]);
         %ylim([-0.1 0.1]);
         
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
         tx= strcat('y= ', '  ', equation, ' (', 'fitted red line', ')');
         text(20, 45, tx, 'FontSize', 10); 
         tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(p_slope(1)),')',...
            ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
         text(20, 40, tx, 'FontSize', 10); 
         tx = strcat('Bias', '(', num2str(p_slope(2)), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
         text(20, 35, tx, 'FontSize', 10); 
         
         %For Reflectance
%          tx= strcat('y= ', '  ', equation, ' (', 'fitted red line', ')');
%          text(245, 0.09, tx, 'FontSize', 10); 
%          tx = strcat('95% Confidence bound for',  {' '}, 'Gain', '(', num2str(p_slope(1)),')',...
%             ':', {' '},'[',  num2str(ci(1,1)), {' '},  {' '}, num2str(ci(2,1)), ']', ',', {' '});
%          text(245, 0.080, tx, 'FontSize', 10); % for Reflectance
%          tx = strcat('Bias', '(', num2str(p_slope(2)), ')', ':',{' '},'[', num2str(ci(1,2)), {' '}, {' '},  num2str(ci(2,2)), ']');
%          text(245, 0.070, tx, 'FontSize', 10); % for Reflectance
         
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