function [Out] = ClusterIntersectHyperion %(search_path)
for cluster = [1,4,13,2:3,5:12,14:19]
    if ~exist(['out',num2str(cluster,'%02d'),'.mat'],'file') % any completed runs exsist 
        display(cluster)
        if ~exist(['out.running_results.C',num2str(cluster,'%02d'),'.mat'],'file') %If not previous partial runs exsit, start from scratch
            clear out
            out.Meta = [];
            out.ZonalBinaryMaskFileName=[];
            save(['out',num2str(cluster,'%02d'),'.mat'],'out')
            out.sensor = 'Hyperion';
            out.sat = 'EO1';
            out.path = 'P181';
            out.row = 'R042';
            out.ClusterNumber=cluster;
            %out.Meta.LPGSVer = [];
            out.ThresholdPixelNumber = 50;
            out.startfile = 2; 
            if ispc
                out.OldPath = addpath('Z:\SpecialNeeds\[Useful Matlab functions]\scan_diamondback');
            elseif isunix
                out.OldPath = addpath('/home/leighl/zdrive/SpecialNeeds/[Useful Matlab functions]/scan_diamondback/');
                %/home/leighl/zdrive/SpecialNeeds/[Useful Matlab functions]/scan_diamondback/
            else
                disp([datstr(now),': error'])
            end
            
            folders=0;
            for paths=160:210
                for rows=36:50
                    if isunix
                        home=getenv('HOME');
                        search_path=fullfile(home,'zdrive','ImageDrive',out.sensor,out.sat,['P',num2str(paths,'%03d')],['R',num2str(rows,'%03d')]);
                        if exist(search_path,'dir')
                            folders=folders+1;
                            out = path_good(out,search_path,folders,paths,rows);
                        end
                    elseif ispc
                        search_path=fullfile('z:','ImageDrive',out.sensor,out.sat,['P',num2str(paths,'%03d')],['R',num2str(rows,'%03d')]);
                        if exist(search_path,'dir')
                            folders=folders+1;
                            out = path_good(out,search_path,folders,paths,rows);
                        end
                    end
                    
                end
            end
        else
            load(['out.running_results.C',num2str(cluster,'%02d'),'.mat'])
            out.startfile = size(out.CountPixel,2)-1; %if a partial run exists start where it ended (minus 1)
            save(['out',num2str(out.ClusterNumber,'%02d'),'.mat'],'out')
        end
        out=Intersect(out);
        save(['out',num2str(out.ClusterNumber,'%02d'),'.mat'],'out')
        %path(out.OldPath)
        Out(cluster)=out;
    end
end
clear global
end

function out=path_good(out,search_path,folders,paths,rows)
    out.path(folders,:) = ['P',num2str(paths,'%03d')];
    out.row(folders,:)  = ['R',num2str(rows, '%03d')];
    current=search_path;
    s = dir(current); 
    s = s(3:end,:);
                    
    [out] = findMetaHyperion(out,s,current);
end

function [out] = findMetaHyperion(out,s,current)
% since hyperion has different name conventions, decided to split the code,
% could maybe of kept is a combined fuction, but I'm sure that would get
% very ugly.  So a lot of repeat of the above with minor differences for
% hyperion
    if isfield(out.Meta, 'LPGSVer')
        count=size(out.Meta.LPGSVer,1);   % LL
    else
        count = 0;
    end
    for i = 1:size(s, 1)
        if (s(i).isdir && regexpi(s(i).name, '^[0-9]{8}$'))
            subcurrent = fullfile(current, s(i).name);
            fprintf('Scanning the folder %s...\n', subcurrent);
            t = dir(subcurrent); t = t(3:end,:);
            fprintf('Looking for the version L1T...\n');
            if (~isL1T(t))
                fprintf('No L1T version...\n');
            %else
                continue;            
            end
            count=count+1;  %LL
            out.date(count,:) = [s(i).name(1:4) '-' s(i).name(5:6) '-' s(i).name(7:8)];
            fprintf('Version L1T found...\n');
            out.level(count,:) = 'L1T'; %LL
            fullpath = fullfile(current, s(i).name, out.level(count,:));
            fprintf('Beginning of the search for the MTL.txt file...\n');
            [u, localpath] = findFileHL1T(fullpath, out.sat);
            fullpath = localpath;
            if (size(u, 1) == 0)
                fprintf('No MTL.txt file found...\n')
                out.meta_file{count} = 'NotAFile';
            error(['no Meta file for L1T image, file = ',s(i).name])
            else
                fprintf('MTL.txt file found...\n');
                out.meta_file{count} = u(1).name;
                meta_fullpath = fullfile(fullpath, str2mat(out.meta_file{count})); %#ok<DSTRMT>
                display([datestr(now),':C',num2str(out.ClusterNumber,'%02d'),': ',meta_fullpath])
                out.Meta = GetMeta(meta_fullpath, out.Meta, out.sat, count);
                out.Meta.Meta_FullPath(count,1)={meta_fullpath};
                out.fullpath{count,1}=fullpath;
            end
            %end
        end
    end
%     folder = fullfile('.', 'output');
%     if (exist(folder, 'dir') == 0)
%         mkdir('output');
%     end
end