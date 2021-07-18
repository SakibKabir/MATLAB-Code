        %%%changes in beta4 
        % 1. added progress display
        clc
        clear all
        
        data_source='highres_kmls';% put the high res kml files here
        data_dest='parent_child_mats';
        files=dir(data_source);
        
        
       %%%create fileby file acess loop
       for i=[4 5]%:length(files)
        
        clearvars -except i data_source data_dest files 
        
        file_path=fullfile(data_source,files(i).name);%%%access the individual files
        output_file=fullfile(data_dest,files(i).name);
        output_file=output_file(:,1:end-4);
        polydata=kml2struct_ndl('S2A_OPER_Grid.kml');
        poly_index=[];

        for n=1:length(polydata)
           poly_index(n)=str2num(polydata(n).Name);%index the good and bad polygon serially 
        end
        
          good_index=find(poly_index);%find the good index ie.,parent
          bad_index=find(~poly_index);%find the bad indexes ie.,child/holes
          parent_counter=0;
              
        for j=good_index %parents
            
                parent_counter=parent_counter+1;
                polygons(parent_counter).parent.lon=polydata(j).Lon(1:end-1);
                polygons(parent_counter).parent.lat=polydata(j).Lat(1:end-1);
                internal_child_counter=0;   
                for k=bad_index %child/holes

                     [in,on] = inpolygon(polydata(k).Lon, polydata(k).Lat ,polydata(j).Lon, polydata(j).Lat);% j is hole and i is parent; check if j is in the i                   
                        if all(in(1:end-1))
                             internal_child_counter= internal_child_counter+1;
                             polygons(parent_counter).child(internal_child_counter).lat=polydata(k).Lat(1:end-1);
                             polygons(parent_counter).child(internal_child_counter).lon=polydata(k).Lon(1:end-1);
                             
                        end
                        %%%diplay the progress
% %                         clc
% %                         disp(['Current progress: ' 'Files:' num2str(i-2) '/' num2str(length(files)-2)...
% %                                           '    Parents:' num2str(round((parent_counter/length(good_index))*100,4)) '% ' '    Current child:' num2str(round((internal_child_counter/length(in))*100,4)) '%' ]);

                        
                end 
               %%%diplay the progress
               clc
               disp(['Current progress: ' 'Files:' num2str(i-2) '/' num2str(length(files)-2)...
                                  '    Parents:' num2str(round((parent_counter/length(good_index))*100,4)) '%']);
              
                
       end
       save(output_file,'polygons');
       disp('done')
       end
       disp('all done')