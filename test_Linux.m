clear all
s = 100;
j=1;
for i = 1:0.1:100
    s1(j) = s+i;
    j = j+1;
end
clear s j i
save('s1')