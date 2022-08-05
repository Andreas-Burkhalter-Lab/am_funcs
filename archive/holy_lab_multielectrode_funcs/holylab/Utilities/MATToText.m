function MATToText(filename)
load(filename)
v = whos;
nvars = size(v);
[varnames{1:nvars}] = deal(v.name);
i = 1;
while ~strcmp('filename',varnames{i})
        i = i+1;
end
keep = setdiff(1:nvars,i);
varnames = varnames(keep);
v = v(keep);
v.class
nvars = nvars-1;
for i = 1:nvars
        temp = eval(varnames{i});
        size(temp);
        dlmwrite([varnames{i},'.txt'],temp);
end
