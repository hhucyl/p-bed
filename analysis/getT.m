function T = getT(yy)
if(numel(yy)>2)
    y = fft(yy);
    y(1) = [];
    n = length(y);
    power = abs(y(1:floor(n/2))).^2;
    maxfreq = 0.5;
    freq = (1:n/2)/(n/2)*maxfreq;
    period = 1./freq;
    [mp,index] = max(power);
    T = period(index);
else
    T =0;
end
end