function [base_L8, date_L8, base_D, date_D] = D_L8_Locations(location)
%%  
if location == 1, Path = 'P039'; Row = 'R037'; elseif location == 2, Path ='P026'; Row ='R029';   
    elseif location == 3, Path ='P021'; Row ='R032'; elseif location == 4, Path = 'P023'; Row ='R034';
    elseif location == 5, Path ='P023'; Row ='R033'; elseif location == 6, Path = 'P029'; Row ='R028';
    elseif location == 7, Path ='P029'; Row ='R029'; elseif location == 8, Path = 'P170'; Row ='R078';
    elseif location == 9, Path ='P170'; Row ='R079'; elseif location == 10,Path = 'P171'; Row ='R078';
    elseif location == 11,Path ='P171'; Row ='R079'; elseif location == 12,Path = 'P174'; Row ='R084';       
    elseif location == 13,Path ='P175'; Row ='R083'; elseif location == 14,Path = 'P072'; Row ='R087';
    elseif location == 15,Path ='P089'; Row ='R083'; elseif location == 16,Path = 'P090'; Row ='R082';
    elseif location == 17,Path ='P091'; Row ='R082'; elseif location == 18,Path = 'P093'; Row ='R086';
    elseif location == 19,Path ='P181'; Row ='R040';
end

if ispc
    % Dove file location and Date Extraction
    base_D = fullfile('Z:\ImageDrive\PlanetLabs\Processed\1047\', Path, Row); 
    date_D = dir(base_D); date_D([1 2])=[]; date_D = date_D.name;

    % Landsat 8 file location and Date Extraction
    base_L8 = fullfile('Z:\ImageDrive\OLI.TIRS\L8\', Path, Row); dates_L8 = dir(base_L8);
    dates_L8([1 2])=[];

    for i = 1: size(dates_L8, 1)
      Date_com = strcmp(dates_L8(i).name, date_D);
      if Date_com == 1
          date_L8 = dates_L8(i).name;
      else   
      end
    end 
elseif isunix
    % Dove file location and Date Extraction
    base_D = fullfile('/home/sakib.kabir/zdrive/ImageDrive/PlanetLabs/Processed/1047/', Path, Row);
    date_D = dir(base_D); date_D([1 2])=[]; date_D = date_D.name;
    
    % Landsat 8 file location and Date Extraction
    base_L8 = fullfile('/home/sakib.kabir/zdrive/ImageDrive/OLI.TIRS/L8/', Path, Row);
    dates_L8 = dir(base_L8); dates_L8([1 2])=[];
   
    for i = 1: size(dates_L8, 1)
      Date_com = strcmp(dates_L8(i).name, date_D);
      if Date_com == 1
          date_L8 = dates_L8(i).name;
      else   
      end
    end 
end
 
 % if you want to write it for linux
%     if ispc
%             
%         elseif isunix
%             base_L8 = '/home/sakib.kabir/zdrive/ImageDrive/OLI.TIRS/L8/P072/R087';
%             dates_L8 = dir(base_L8);
%         end
% 
%         % Dove file location
%         if ispc
%             base_D ='Z:\ImageDrive\PlanetLabs\Processed\1047\P072\R087'; 
%             date_D = dir(base_D); date_D([1 2])=[]; 
%         elseif isunix
%             base_D ='/home/sakib.kabir/zdrive/ImageDrive/PlanetLabs/Processed/1047/P072/R087'; 
%             date_D = dir(base_D);
%         end  
%         
%      elseif location == 2
% 
% end
% end


