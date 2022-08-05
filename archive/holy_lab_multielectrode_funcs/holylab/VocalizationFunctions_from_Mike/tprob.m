table = [codes(1:end-1)',codes(2:end)'];
intraboutidx = find(p < cutoff);
newtable = table(intraboutidx,1:2);
nbins = {1:1:91,1:1:91};
Z = hist3(newtable,nbins);
Z = Z./sum(sum(Z));