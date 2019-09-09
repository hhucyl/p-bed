clear
clc
prefix_name = {'/home/user/p-bed/'};
name = strcat(prefix_name,'time.txt');
N = 87600;
fid = fopen(char(name));
t = [];
for i=1:N*2
    if(mod(i,2)==0)
    else
        ii = (i+1)/2;
        a{ii,1} = str2num(fgets(fid));
        ii
        temp = a{ii,1};
        t = [t;temp(2:end)'];
    end
end
hist(t)