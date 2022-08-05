function sniptime(idx,snips,twhis)

x = snips{idx};
x = abs(x);
x = sum(x,1);
k = find(x);

t = linspace(twhis(1,idx),twhis(2,idx),size(x,2));
tk = t(k);
y = log10(x(k));

plot(tk,y);







end
