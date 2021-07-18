%simply generates a fractional day of week output 
% input is a datenum
% Output is Monday = 0 , Sunday = 6
% Decimal compponet is hour fraction of the day. 
function DF=DayFrac(datenumber)
[~,~,~,h,mn,~]=datevec(datenumber);
DF=mod(weekday(datenumber+5),7)+((h+mn/60)/24);