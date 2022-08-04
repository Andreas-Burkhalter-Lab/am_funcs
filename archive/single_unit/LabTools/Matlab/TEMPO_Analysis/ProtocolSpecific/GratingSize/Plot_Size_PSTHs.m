lines = {'b' 'r' 'g'};
lines2 = {'b-.' 'r-.' 'g-.'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non-Normalized Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
hold on
plot(0,0,'b'); plot(0,0,'b-.');
plot(0,0,'r'); plot(0,0,'r-.');
plot(0,0,'g'); plot(0,0,'g-.');
legend('V1,StimType1','V1,StimType2','LM,StimType1','LM,StimType2','AL,StimType1','AL,StimType2')
title_str = ['PSTH at Optimized Size']
title(title_str),xlabel('Time (s)'), ylabel('Spikes')

for n = 1:3;
    if n == 1;
        OptimizedSize_PSTH = load('V1OptimizedSize_PSTH.dat');
        OptimizedSize_PSTH2 = load('V1OptimizedSize_PSTH2.dat');
    elseif n == 2;
        OptimizedSize_PSTH = load('LMOptimizedSize_PSTH.dat');
        OptimizedSize_PSTH2 = load('LMOptimizedSize_PSTH2.dat');
    elseif n == 3;
        OptimizedSize_PSTH = load('ALOptimizedSize_PSTH.dat');
        OptimizedSize_PSTH2 = load('ALOptimizedSize_PSTH2.dat');
    end

    [rows, nbins] = size(OptimizedSize_PSTH);
    bins = OptimizedSize_PSTH(1,:);
    total_cells = rows-1
    hits = [];
    for j = 2:nbins;
        hits(j) = sum(OptimizedSize_PSTH(:,j));
    end
    hits = hits/total_cells;
    hold on
    plot(bins,hits,lines{n})


    [rows, nbins] = size(OptimizedSize_PSTH2);
    bins = OptimizedSize_PSTH2(1,:);
    total_cells = rows-1
    hits = [];
    for j = 2:nbins;
        hits(j) = sum(OptimizedSize_PSTH2(:,j));
    end
    hits = hits/total_cells;
    hold on
    plot(bins,hits,lines2{n})

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalized Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
hold on
plot(0,0,'b'); plot(0,0,'b-.');
plot(0,0,'r'); plot(0,0,'r-.');
plot(0,0,'g'); plot(0,0,'g-.');
legend('V1,StimType1','V1,StimType2','LM,StimType1','LM,StimType2','AL,StimType1','AL,StimType2')
title_str = ['PSTH at Optimized Size (Cell Responses Normalized)']
title(title_str),xlabel('Time (s)'), ylabel('Spikes')

for k = 1:3;

    if k == 1;
        OptimizedSize_PSTH_Norm = load('V1OptimizedSize_PSTH_Norm.dat');
        OptimizedSize_PSTH_Norm2 = load('V1OptimizedSize_PSTH_Norm2.dat');
    elseif k == 2;
        OptimizedSize_PSTH_Norm = load('LMOptimizedSize_PSTH_Norm.dat');
        OptimizedSize_PSTH_Norm2 = load('LMOptimizedSize_PSTH_Norm2.dat');
    elseif k == 3;
        OptimizedSize_PSTH_Norm = load('ALOptimizedSize_PSTH_Norm.dat');
        OptimizedSize_PSTH_Norm2 = load('ALOptimizedSize_PSTH_Norm2.dat');
    end

    [rows, nbins] = size(OptimizedSize_PSTH_Norm);
    bins = OptimizedSize_PSTH_Norm(1,:);
    total_cells = rows-1
    hits = [];
    for j = 2:nbins;
        hits(j) = sum(OptimizedSize_PSTH_Norm(:,j));
    end
    hits = hits/max(hits);
    hold on
    plot(bins,hits,lines{k})


    [rows, nbins] = size(OptimizedSize_PSTH_Norm2);
    bins = OptimizedSize_PSTH_Norm2(1,:);
    total_cells = rows-1
    hits = [];
    for j = 2:nbins;
        hits(j) = sum(OptimizedSize_PSTH_Norm2(:,j));
    end
    hits = hits/max(hits);
    hold on
    plot(bins,hits,lines2{k})
end


