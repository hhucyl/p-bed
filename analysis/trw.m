clear
clc
prefix = {'/home/user/p-bed/'};
N = 1e5;
D = 0.01;
t = 1e4;
name = strcat(prefix,'test_rw1_',num2str(99,'%04d'),'.h5');
con = double(h5read(char(name),char('/Con')));
nx = h5read(char(name),char('/Nx'));
ny = h5read(char(name),char('/Ny'));
con = reshape(con,[nx,ny]);
plot(con(50,:)./N,'*')
hold on
xx = 1:double(nx);
c = 1/(4*D*pi*t).*exp(-abs(xx-50).^2./(4*D*t));
plot(xx,c,'r.')