clear
clc
prefix_name = {'/home/user/p-bed/'};
num = [0:999];
ii = 1;
for i=1:numel(num)
    
    name = strcat(prefix_name{ii},'test_pbed_r1_',num2str(num(i),'%04d'),'.h5');
    nx = h5read(char(name),char('/Nx'));
    ny = h5read(char(name),char('/Ny'));
    v(:,i) = h5read(char(name),char('/Velocity_0'));
    i
end
va = mean(v,2);
for i = 1:numel(num)
    Vhas(:,i) = v(:,i)-va;
    uhx(:,:,i) = reshape(Vhas(1:3:end-2,i),[nx,ny]);
    vhy(:,:,i) = reshape(Vhas(2:3:end-1,i),[nx,ny]);
    uvh(:,:,i) = uhx(:,:,i).*vhy(:,:,i);
    i
    
end
