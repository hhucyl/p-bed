P = [0.002747,0.006077,0.007429,0.008912];
Y = [43,65,86,108];
N = [20,30,40,50]
V = 220.*(Y-1)-N.*pi*8*8;
Dm = 2.3e-6;
y= pi.*(0.5.*V./220.*P).^2./Dm
plot(y,Y./160-108/160,'*')
ylabel('Depth')
xlabel('\itD_{\rmeff}/\itD_{\rmm}')