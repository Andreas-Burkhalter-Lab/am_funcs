function snipplot(idx,snips,twhis)

x = snips{idx};
[i,j,k] = find(x);
k = abs(k);
k = log(k);
k = k./min(k);
y = sparse(i,j,k,size(x,1),size(x,2));

t = linspace(twhis(1,idx),twhis(2,idx),size(y,2));
f = linspace(0,125,257);

h = pcolor(t,f,full(y));
set(h, 'EdgeColor', 'none');
colormap(1-gray);

end
