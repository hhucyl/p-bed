clear
clc
name = {'/home/pzhang/chen/p-bed/analysis/M.h5'};
M = h5read(char(name),char('/M'));
m = M(5,:)./M(5,1);
t = sqrt(0:999);
kkk = 1:20*20;
plot(t,m,'*')
hold on
tt = t(kkk);
mm = m(kkk);
plot(t(kkk),m(kkk),'r*')