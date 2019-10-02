clc
m = size(pdy);
TT = [];
YY = [];
for i=1:m(1)
for j=1:m(2)
if(numel(pdy{i,j})>0)
yy = abs(diff(pdy{i,j}));
% nx = 224;
% yy(yy>100) = abs(yy(yy>100)-nx);

% subplot(3,1,1)
% plot(yy,'*')
y = fft(yy);
y(1) = [];
n = length(y);
power = abs(y(1:floor(n/2))).^2;
maxfreq = 0.5;
freq = (1:n/2)/(n/2)*maxfreq;
period = 1./freq;
% kkk = find(period<180);
[mp,index] = max(power);
T = period(index);
TT = [TT;T];
YY = [YY;pdy{i,j}'];

% subplot(312)
% plot(freq,power)
% subplot(313)
% plot(period,power)
% title(i)
% drawnow
clf
i
end
end
end
hist(TT(TT>0),[0:50:700 800:500:max(TT)])