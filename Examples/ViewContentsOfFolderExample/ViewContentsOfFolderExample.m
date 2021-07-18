%% View Contents of Folder
% List the contents of a folder.
%
% Create a folder, |myfolder|, that contains the files |myfile1.m|, 
% |myfile2.m|, and |myfile3.m|.

mkdir myfolder
movefile myfile1.m myfolder
movefile myfile2.m myfolder
movefile myfile3.m myfolder

%% 
% List the files in |myfolder|.
dir myfolder
