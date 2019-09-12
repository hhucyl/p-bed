nn = 1;
TT = [];
m = size(pdy);
for i=1:m(1)
for j=1:m(2)
if(numel(pdy{i,j})>0)
    TT(nn,1) = getT(pdy{i,j});
    nn = nn+1;
    i
end
end
end
t = [0,50,100,200,300,400:200:1e3];
hist(TT,t)