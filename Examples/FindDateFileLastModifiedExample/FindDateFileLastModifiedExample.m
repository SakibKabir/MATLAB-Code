%% Find Date File Last Modified
% Get the serial date number for the date and time a file was last modified.
%
% Use the |datenum| field of the structure returned by the |dir| command.
% Do not use the |datenum| function to convert the |date| field of the 
% structure to a number. The results of the |datenum| function vary depending 
% on the locale. Instead, use the |datenum| field.
MyFileInfo = dir('myfile1.m');
FileDate = MyFileInfo.datenum