clear
clc
prefix_name = {'/home/user/p-bed/'};
prefix_name = {'/media/user/PZ_Q/p-bed/1/'};
prefix_name = {'/media/user/PZ_Q/p-bed/1e4/'};
fname = {'test_pbed_t1_','test_pbed_t2_'};

color = {'r','b'};
dt = [1e3,1e3];
for ii = 2:2
pp = h5read(char(strcat(prefix_name,fname{ii},num2str(0,'%04d'),'.h5')),'/RWPposition');
NN = numel(pp)/3;
N = 87600;
ppy = 130;
num = [0:2999];
t = zeros(N,numel(num));
Px = zeros(N,numel(num));
Py = zeros(N,numel(num));
for i=1:numel(num)
    name = strcat(prefix_name,fname{ii},num2str(num(i),'%04d'),'.h5');
    p = h5read(char(name),char('/RWPposition'));
    px = p(1:3:end-2);
    py = p(2:3:end-1);
    Px(:,i) = px;
    Py(:,i) = py;
    kkk = find(py<ppy);
    t(kkk,i) = 1;
    i
end
T = [];
for j = 1:N
    k1 = find(t(j,:)==1);
    temp = mat2cell(k1,1,diff([0,find(diff(k1)~=1),length(k1)]));
    ttemp= [];
    tktemp = [];
    for jj=1:numel(temp)
        if(numel(temp{jj}).*dt(ii)==0)
            tttemp = [];
            tkttemp = [];
        else
            tttemp = numel(temp{jj}).*dt(ii);
            tkttemp = min(temp{jj});
        end
        ttemp = [ttemp,tttemp];
        tktemp = [tktemp,tkttemp];
    end
    tt{j,1} = ttemp;
    tk{j,1} = tktemp;
    T = [T;tt{j,1}'];
    TT{ii} = T;
    j
end
ttt = 1e3:1e3:max(T);
y = hist(T,ttt)./numel(T);
figure(1)
loglog(ttt,y,char(strcat(color{ii},'*')))
hold on
figure(2)
h = y./numel(T);
cy = cumsum(h)/sum(h);
plot(ttt,cy,char(strcat(color{ii},'*')))
hold on
set(gca,'xscale','log')
drawnow
end

% index = [];
% for i=1:N
%     kk1 = find(tt{i}(1:end)>5e3);
%     kk2 = find(tt{i}(1:end)>1e4);
%     if(numel(kk1)>0 && numel(kk2)==0)%(numel(kk)==0)
%         index = [index;i];
%     end
%     
% end
% 
% num = [0:999];
% figure(3)
% name = strcat(prefix_name,'test_pbed_t1_',num2str(0,'%04d'),'.h5');
% nx = double(h5read(char(name),char('/Nx')));
% ny = double(h5read(char(name),char('/Ny')));
% for i=1:numel(num)
%     name = strcat(prefix_name,'test_pbed_t1_',num2str(num(i),'%04d'),'.h5');
%     p = h5read(char(name),char('/RWPposition'));
%     px = p(1:3:end-2);
%     py = p(2:3:end-1);
%     plot(px(index),py(index),'b*')
%     hold on
%     plot([0 nx], [109 109],'r')
%     axis([0 nx 0 ny])
%     drawnow
%     clf
% end

