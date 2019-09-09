clear
clc
prefix_name = {'/home/user/p-bed/'};
num = [0:999];
cx = 112;
cy = [25 45 65 85];
color = {'r','b','g','k'};
name = strcat(prefix_name,'test_pbed_r_',num2str(0,'%04d'),'.h5');
nx = double(h5read(char(name),char('/Nx')));
ny = double(h5read(char(name),char('/Ny')));
p = h5read(char(name),char('/RWPposition'));
px = p(1:3:end-2);
py = p(2:3:end-1);
for i=1:numel(cy)
    r = sqrt((px-cx).^2+(py-cy(i)).^2);
    kkk{i} = find(r<=5);
end
for i=1:numel(num)
    name = strcat(prefix_name,'test_pbed_r_',num2str(num(i),'%04d'),'.h5');
    p = h5read(char(name),char('/RWPposition'));
    px = p(1:3:end-2);
    py = p(2:3:end-1);
    %hist(py)
    for j=1:numel(cy)
        ppx = px(kkk{j});
        ppy = py(kkk{j});
        M = numel(ppx);
        Xc = sum(ppx)/M;
        Yc = sum(ppy)/M;
        ax1 = subplot(3,1,1)
        plot(px(kkk{j}),py(kkk{j}),strcat(char(color(j)),'.'))
        hold on
        axis([0 nx 0 ny])
        subplot(3,1,2)
        plot(Xc,Yc,strcat(char(color(j)),'*'))
        axis([0 nx 0 ny])
        hold on
        subplot(313)
        oyy = sum(ppy.^2)/M-Yc^2;
        oxy = sum(ppy.*ppx)/M-Yc*Xc;
        plot(i,oyy,strcat(char(color(j)),'*'))
        hold on
    end
    
    drawnow
    title(i)
    %saveas(gcf,char(strcat(prefix_name,num2str(i),'.png')));
    %pause(1)
    cla(ax1)
end