% Decide how many groups there are:

delta = f3(2:end) - f3(1:end-1);
k = find(delta >= 2);
groups = length(k) + 1;
f4 = zeros(1,groups);

tbin = 1:1:max(t);
ft = cell(1,tbin);
sizes = zeros(1,tbin);
for i = 1:length(tbin)
    tp = t(t==tbin(i));
    fp = f2(t==tbin(i));
    if length(fp) == 1
        ft{i} = fp;
    elseif length(fp) > 1
        delta = abs(fp(2:end) - fp(1:end-1));
        k = find(delta >= 10);
        groups = [1,((k)'+1)];
        f5 = fp(groups); % initialize frequencies to test
        for j = 1:length(f5)
            a = find(abs(fp - f5(j)) <= 10);
            f5(j) = mean(fp(a));
        end
        ft{i} = f5;
    end
    sizes(i) = length(ft{i});
end

sizebins = 1:1:max(sizes);

for i = 1:length(sizebins);
    h{i} = find(sizes == sizebins(i));
    x = ft(h{i});
    fh{i} = cat(2,x{:});
end

% Edit this part. Harms should take a harmonic only if there are more that
% 5 msec of continuous overlap
for i = 2:length(fh) % All tbins with 2 or more frequencies
    b = fh{i};
    c = h{i};
    if size(b,2) >= 5 % At least 5 msec of harmonic overlaps
        harms{i-1} = zeros(1,i);
        for j = 1:length(harms{i-1})
            harms{i-1}(j) = mean(b(j,:));
        end
    end
end
    
        