%need to run ReadTime first
clc
tt = min(t):50:max(t);
y = hist(t,tt);
figure(1)
loglog(tt,y./N,'*')
figure(2)
ttt = cumsum(y./N);
cdfplot(t)
index = [];
for i=1:N
    kk = find(a{i}(2:end)>1e4);
    if(numel(kk)==0)%(numel(kk)==0)
        index = [index;a{i}(1)+1];
    end
    
end

num = [0:999];
figure(3)
name = strcat(prefix_name,'test_pbed_t1_',num2str(0,'%04d'),'.h5');
nx = double(h5read(char(name),char('/Nx')));
ny = double(h5read(char(name),char('/Ny')));
for i=1:numel(num)
    name = strcat(prefix_name,'test_pbed_t1_',num2str(num(i),'%04d'),'.h5');
    p = h5read(char(name),char('/RWPposition'));
    px = p(1:3:end-2);
    py = p(2:3:end-1);
    plot(px(index(:)),py(index(:)),'b*')
    hold on
    plot([0 nx], [109 109],'r')
    axis([0 nx 0 ny])
    drawnow
    clf
end