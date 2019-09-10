clear
clc
prefix = {'/home/user/p-bed/'};
name = strcat(char(prefix),'test_pbed1_0001.h5');
pos = h5read(char(name),'/Pposition');
N = numel(pos)/6;
px = pos(1:3:3*N-2);
py = pos(2:3:3*N-1);
plot(px,py,'*')
outname = strcat(char(prefix),'circle_p.txt');
fid = fopen(char(outname),'w');
fprintf(fid,'%f\n',0.5);
fprintf(fid,'%d\n',N);
for i = 1:N
    fprintf(fid,'%f %f %f\n',[px(i),py(i) 10.2]);
end
fclose(fid)