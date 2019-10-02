clc
%t_r
index = [];
ik = [];
for i=1:N
    kk1 = find(tt{i}(1:end)>2.2655e4);
    kk2 = find(tt{i}(1:end)<=1e7);
    kk = intersect(kk1,kk2);
    if(numel(kk)>0)%(numel(kk)==0)
        index = [index;i];
        ii = numel(index);
        ik{ii,1}= tk{i}(kk);%get in time
        ik{ii,2}= tk{i}(kk)+tt{i}(kk)./1e3-1;%get out time
    end
    i
    
end
PDY = [];
PDX = [];
DPDY = [];
DPDX = [];
pdx=[];
pdy=[];
for i=1:numel(index)
    tau = [];
    tau(:,1) = ik{i,1}';
    tau(:,2) = ik{i,2}';
    for j=1:(numel(tau)/2)
        pdx{i,j} = Px(index(i),tau(j,1):tau(j,2));
        pdy{i,j} = Py(index(i),tau(j,1):tau(j,2));
        %PDY = [PDY;pdy{i,j}'];
        %DPDY = [DPDY;diff(pdy{i,j}')];
        %PDX = [PDX;pdx{i,j}'];
        %DPDX = [DPDX;diff(pdx{i,j}')];
%         plot(pdx{i,j},pdy{i,j},'b*')
%         hold on
%         axis([0 220 0 170])
%         
%         drawnow
    end
%     title(i)
%     clf
    
    i
end

% num = [0:999];
% figure(3)
% name = strcat(prefix_name,'test_pbed_t1_',num2str(0,'%04d'),'.h5');
% nx = double(h5read(char(name),char('/Nx')));
% ny = double(h5read(char(name),char('/Ny')));
% ppy = 109;
% pdd = [];
% for i=1:numel(num)
%     name = strcat(prefix_name,'test_pbed_t1_',num2str(num(i),'%04d'),'.h5');
%     p = h5read(char(name),char('/RWPposition'));
%     px = p(1:3:end-2);
%     py = p(2:3:end-1);
%     for j=1:numel(index)
%         if
%         px(index(j))
%     end
%     plot(px(index(1)),py(index(1)),'b*')
% %     pd = (py(index)-ppy);
% %     pdd = [pdd;pd(pd<0)];
%     hold on
%     plot([0 nx], [109 109],'r')
%     axis([0 nx 0 ny])
%     title(i)
%     drawnow
%     clf
% end