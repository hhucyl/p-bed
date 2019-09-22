clear
clc
num = [594:673];
R = 8;
for i = 1:numel(num)
    name = strcat('test_pbed_r1_',num2str(num(i),'%04d'),'.h5');
    p = h5read(char(name),char('/RWPposition'));
    pb = h5read(char(name),char('/RWPpositionb'));
    px = p(1:3:end-2);
    py = p(2:3:end-1);
    pxb = pb(1:3:end-2);
    pyb = pb(2:3:end-1);
    pos = h5read(char(name),char('/Pposition'));
    
    Np = numel(pos)/6;
    posx = pos(1:3:end-2);
    posy = pos(2:3:end-1);
    kkk = [];
    for j = 1:Np
        l = sqrt((px-posx(j)).^2+(py-posy(j)).^2);
        kkk = [kkk,find(l<R)];
    end
    count(i,1) = numel(kkk);
    RR = zeros(Np,1)+R;
    viscircles([posx(1:Np),posy(1:Np)],RR);
    hold on
%     pp(i,:) =[px(23241),py(23241)];
%     ppb(i,:) =[pxb(23241),pyb(23241)]; 
%     plot(pp(:,1),pp(:,2),'-k*')
%     plot(ppb(:,1),ppb(:,2),'b*')
    plot(px(kkk),py(kkk),'k*')
    axis equal
    title(num(i))
    drawnow
%     clf
end