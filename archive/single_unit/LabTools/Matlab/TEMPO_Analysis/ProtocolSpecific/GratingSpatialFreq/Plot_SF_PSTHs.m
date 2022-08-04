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
title_str = ['PSTH at Optimized Spatial Freq']
title(title_str),xlabel('Time (s)'), ylabel('Spikes')

for n = 1:3;
    if n == 1;
        OptimizedSF_PSTH = load('V1OptimizedSF_PSTH.dat');
        OptimizedSF_PSTH2 = load('V1OptimizedSF_PSTH2.dat');
    elseif n == 2;
        OptimizedSF_PSTH = load('LMOptimizedSF_PSTH.dat');
        OptimizedSF_PSTH2 = load('LMOptimizedSF_PSTH2.dat');
    elseif n == 3;
        OptimizedSF_PSTH = load('ALOptimizedSF_PSTH.dat');
        OptimizedSF_PSTH2 = load('ALOptimizedSF_PSTH2.dat');
    end

    [rows, nbins] = size(OptimizedSF_PSTH);
    bins = OptimizedSF_PSTH(1,:);
    total_cells = rows-1
    hits = [];
    for j = 2:nbins;
        hits(j) = sum(OptimizedSF_PSTH(:,j));
    end
    hits = hits/total_cells;
    hold on
    plot(bins,hits,lines{n})


    [rows, nbins] = size(OptimizedSF_PSTH2);
    bins = OptimizedSF_PSTH2(1,:);
    total_cells = rows-1
    hits = [];
    for j = 2:nbins;
        hits(j) = sum(OptimizedSF_PSTH2(:,j));
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
title_str = ['PSTH at Optimized Spatial Freq (Cell Responses Normalized)']
title(title_str),xlabel('Time (s)'), ylabel('Spikes')

for k = 1:3;

    if k == 1;
        OptimizedSF_PSTH_Norm = load('V1OptimizedSF_PSTH_Norm.dat');
        OptimizedSF_PSTH_Norm2 = load('V1OptimizedSF_PSTH_Norm2.dat');
    elseif k == 2;
        OptimizedSF_PSTH_Norm = load('LMOptimizedSF_PSTH_Norm.dat');
        OptimizedSF_PSTH_Norm2 = load('LMOptimizedSF_PSTH_Norm2.dat');
    elseif k == 3;
        OptimizedSF_PSTH_Norm = load('ALOptimizedSF_PSTH_Norm.dat');
        OptimizedSF_PSTH_Norm2 = load('ALOptimizedSF_PSTH_Norm2.dat');
    end

    [rows, nbins] = size(OptimizedSF_PSTH_Norm);
    bins = OptimizedSF_PSTH_Norm(1,:);
    total_cells = rows-1
    hits = [];
    for j = 2:nbins;
        hits(j) = sum(OptimizedSF_PSTH_Norm(:,j));
    end
    hits = hits/max(hits);
    hold on
    plot(bins,hits,lines{k})


    [rows, nbins] = size(OptimizedSF_PSTH_Norm2);
    bins = OptimizedSF_PSTH_Norm2(1,:);
    total_cells = rows-1
    hits = [];
    for j = 2:nbins;
        hits(j) = sum(OptimizedSF_PSTH_Norm2(:,j));
    end
    hits = hits/max(hits);
    hold on
    plot(bins,hits,lines2{k})
end


