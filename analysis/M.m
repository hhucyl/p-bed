clear
clc
name = {'/home/user/p-bed/analysis/M.h5'};
M = h5read(char(name),char('/M'));
m = M(4,:)./M(4,1);
t = sqrt(0:999);
kkk = 1:20*20;
plot(t,m,'*')
hold on
tt = t(kkk);
mm = m(kkk);
plot(t(kkk),m(kkk),'r*')