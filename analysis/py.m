clear
clc
prefix_name = {'/home/user/p-bed/'};
fname = {'test_pbed_t1_','test_pbed_t2_'};
color = {'r','b'};
num = [0:999];
for ii = 1:2
figure(3)
name = strcat(prefix_name,fname{ii},num2str(0,'%04d'),'.h5');
nx = double(h5read(char(name),char('/Nx')));
ny = double(h5read(char(name),char('/Ny')));
ppy = 109;
pdd = [];
for i=1:numel(num)
    name = strcat(prefix_name,fname{ii},num2str(num(i),'%04d'),'.h5');
    p = h5read(char(name),char('/RWPposition'));
    px = p(1:3:end-2);
    py = p(2:3:end-1);
    dpy = py-ppy;
    ddpy = dpy(dpy<0);
    plot(i,(median(ddpy))./16,char(strcat(char(color{ii}),'*')))
    hold on
    drawnow
end
end