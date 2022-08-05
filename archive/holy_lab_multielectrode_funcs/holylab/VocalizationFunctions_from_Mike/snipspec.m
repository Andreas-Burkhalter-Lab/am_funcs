function snipspec(idx,snips)

x = snips{idx};
x = abs(x);
x = sum(x,2);
k = find(x);

f = linspace(0,125,257);
fk = f(k);
y = log10(x(k));

plot(fk,y);







end
