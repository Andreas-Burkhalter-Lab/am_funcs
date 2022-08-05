function[ujumps,djumps,alljumps] = jumpSizes(pf)

    jumps = [500,1000,5000,10000,15000,20000,25000,30000,35000,40000,45000,50000];
    x = diff(pf(pf>0));
    
    ujumps = zeros(1,length(jumps));
    djumps = zeros(1,length(jumps));
    alljumps = zeros(1,length(jumps));
    for j = 1:length(jumps)
        ujumps(j) = numel(find(x >= jumps(j)));
        djumps(j) = numel(find(x <= -jumps(j)));
        alljumps(j) = numel(find(abs(x) >= jumps(j)));
    end

end