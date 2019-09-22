clear
clc
prefix = {'/home/user/p-bed/'};
num = [0:689];
name = strcat(prefix,'test_pbed_r1_',num2str(0,'%04d'),'.h5');
nx = double(h5read(char(name),'/Nx'));
ny = double(h5read(char(name),'/Ny'));
p = h5read(char(name),char('/RWPposition'));
py = p(2:3:end-1);
k1 = [43,65,86,108];
k1 = [66,88,110,131];
color = {'r','k','c','b'};
ppy = 109;
ppy = 131;
for i = 1:numel(k1)
    kkk{i,1} = find(py<=k1(i));
    
    M0(1,i) = numel(kkk{i});
end
for i=1:numel(num)
    name = strcat(prefix,'test_pbed_r1_',num2str(num(i),'%04d'),'.h5');
    p = h5read(char(name),char('/RWPposition'));
    px = p(1:3:end-2);
    py = p(2:3:end-1);
    %plot(px(kkk{4}),py(kkk{4}),'*')
    %drawnow
    for j=1:numel(k1)
        M(i,j) = numel(find(py(kkk{j,1})<=k1(j)));
        subplot(1,2,1)
        plot(i,M(i,j)./M0(j),char(strcat(color(j),'*')))
        hold on
        subplot(1,2,2)
        if(j==1)
            kkk1 = intersect(find(py<=k1(j)),find(py>=1));
        else
            kkk1 = intersect(find(py<=k1(j)),find(py>=k1(j-1)));
        end
        plot(i,numel(kkk1),char(strcat(color(j),'*')))
        hold on
        
        drawnow
    end
    if(i==1)
        legend('1','2','3','4','location','eastoutside')
    end
end
tt = sqrt(num);