clear
clc
prefix_name = {'/media/user/PZ_Q/p-bed/1/','/home/user/p-bed/analysis/'};
prefix_name = {'/home/user/p-bed/analysis/'};
num = [0:999];color = {'r','b'};
ppy=109;
ppy = 130;
for ii = 1:1
for i=1:numel(num)
    
    name = strcat(prefix_name{ii},'turbulence_',num2str(num(i),'%04d'),'.h5');
    nx = h5read(char(name),char('/Nx'));
    ny = h5read(char(name),char('/Ny'));
    ga = h5read(char(name),char('/Gamma'));
    vel = h5read(char(name),char('/Vhas'));
    vx = reshape(vel(1:3:end-2),[nx,ny]);
    vy = reshape(vel(2:3:end-1),[nx,ny]);
    uhas(:,:,i) = vx;
    vhas(:,:,i) = vy;
    ga(ga>1) = 1.0;
    ga = reshape(ga,[nx,ny]);
    uv(:,:,i) = vx.*vy.*(1-ga);
    uhas(:,:,i) = uhas(:,:,i).*(1-ga);
    vhas(:,:,i) = vhas(:,:,i).*(1-ga);
    ubed(i,1) = sum(sum(uhas(:,1:ppy-1,i).^2));
    vbed(i,1) = sum(sum(vhas(:,1:ppy-1,i).^2));
    
end
vel = h5read(char(name),char('/Vave'));
ua = reshape(vel(1:3:end-2),[nx,ny]);
va = reshape(vel(2:3:end-1),[nx,ny]);
uaa = mean(ua,1);
vaa = mean(va,1);
for i= 1:numel(uaa)
    us(:,i) = ua(:,i) - uaa(i);
    vs(:,i) = va(:,i) - vaa(i);
end


uvs = us.*vs;
figure(1)
plot(mean(uvs,1))
figure(2)
uva = mean(uv,3);
pcolor(uva)
figure(3) %u* Rek
uvhasy = mean(uva,1);
hold on
plot(uvhasy,char(color{ii}))
uu = sqrt(-min(uvhasy))
gga = mean(ga,1);
py = 130;

fi = (py*double(nx)-60*pi*8^2)/(py*double(nx));
k = 5.6e-3*fi^3/(1-fi)^2*16*16;
Rek = uu*sqrt(k)/0.0005
figure(4)
subplot(311)
kbed = 0.5.*(ubed+vbed);
plot(kbed,'*')
y = fft(kbed);
y(1) = [];
n = length(y);
power = abs(y(1:floor(n/2))).^2;
maxfreq = 0.5;
freq = (1:n/2)/(n/2)*maxfreq;
period = 1./freq;
subplot(312)
plot(freq,power)
subplot(313)
plot(period,power)
[mp,index] = max(power);
T = getT(kbed)




end

