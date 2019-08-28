clear
clc
prefix_name = {'/home/user/p-bed/'};
num = [0:999];
for i=1:numel(num)
    
    name = strcat(prefix_name,'analysis/','turbulence_',num2str(num(i),'%04d'),'.h5');
    nx = h5read(char(name),char('/Nx'));
    ny = h5read(char(name),char('/Ny'));
    vel = h5read(char(name),char('/Vhas'));
    vx = reshape(vel(1:3:end-2),[nx,ny]);
    vy = reshape(vel(2:3:end-1),[nx,ny]);
    uv(:,:,i) = vx.*vy;
end
uva = mean(uv,3);

uvhasy = mean(uva,1);
plot(uvhasy)
