function [dec_yr]=str2dec_yr(str)

datenumber=datenum(strcat(str(1:4),'-',str(5:6),'-',str(7:8)));
[~,frac]=date2doy(datenumber);
dec_yr=str2num(str(1:4))+frac;

end