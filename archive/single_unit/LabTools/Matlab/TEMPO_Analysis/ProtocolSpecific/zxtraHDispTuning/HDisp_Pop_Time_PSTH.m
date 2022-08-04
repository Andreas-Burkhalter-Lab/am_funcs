function HDisp_Pop_Time_PSTH(FILE)

load(FILE);

figure;
plot(time(1,:), nanmean(mov_DDI), 'r-')
hold on;
plot(time(1,:), nanmean(stat_DDI), 'r--')
plot(time(1,:), nanmean(mov_resp), 'k-')
plot(time(1,:), nanmean(stat_resp), 'k--')    

%Output = [time(1,:)' nanmean(stat_resp)' nanmean(stat_DDI)' nanmean(mov_resp)' nanmean(mov_DDI)' ];
outfile = [FILE '.out']
%save(outfile, 'Output', '-ascii')

return;

