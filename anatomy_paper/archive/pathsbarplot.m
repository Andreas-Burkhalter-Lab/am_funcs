npaths = 9;

amygrows = [2 3 4];
 lpporrows = [11 12];
 dlgnporrows = [19 20];
 lplmrows = [21];
 dlgnlmrows = [23 24];
 lplirows = [22];
 dlgnlirows = [25 26];
lpv1rows = [6 7 8];
dlgnv1rows = [13 14 15 16 18];

v1toporrows = [27];


amygvals = allpaths.patchInterpatchRatio(amygrows)';
lpporvals = allpaths.patchInterpatchRatio(lpporrows)';
dlgnporvals = allpaths.patchInterpatchRatio(dlgnporrows)';
lplmvals = allpaths.patchInterpatchRatio(lplmrows)';
dlgnlmvals = allpaths.patchInterpatchRatio(dlgnlmrows)';
lplivals = allpaths.patchInterpatchRatio(lplirows)';
dlgnlivals = allpaths.patchInterpatchRatio(dlgnlirows)';
lpv1vals = allpaths.patchInterpatchRatio(lpv1rows)';
dlgnv1vals = allpaths.patchInterpatchRatio(dlgnv1rows)';

v1toporvals = allpaths.patchInterpatchRatio(v1toporrows)';


%standard error
amygster = std(amygvals) / sqrt(length(amygvals));
lpporster = std(lpporvals) / sqrt(length(lpporvals));
dlgnporster = std(dlgnporvals) / sqrt(length(dlgnporvals));
lplmster = std(lplmvals) / sqrt(length(lplmvals));
dlgnlmster = std(dlgnlmvals) / sqrt(length(dlgnlmvals));
lplister = std(lplivals) / sqrt(length(lplivals));
dlgnlister = std(dlgnlivals) / sqrt(length(dlgnlivals));
lpv1ster = std(lpv1vals) / sqrt(length(lpv1vals));
dlgnv1ster = std(dlgnv1vals) / sqrt(length(dlgnv1vals));

amygmean = mean(amygvals);
lppormean = mean(lpporvals);
dlgnpormean = mean(dlgnporvals);
lplmmean = mean(lplmvals);
dlgnlmmean = mean(dlgnlmvals);
lplimean = mean(lplivals);
dlgnlimean = mean(dlgnlivals);
lpv1mean = mean(lpv1vals);
dlgnv1mean = mean(dlgnv1vals);

v1topormean = mean(v1toporvals);

% % % % figure
% % % % % % scatter(1:npaths,[amygmean,lppormean,dlgnpormean,lpv1mean,dlgnv1mean,v1topormean])
% % % % % % scatter(1:npaths,[amygmean,lppormean,dlgnpormean,lpv1mean,dlgnv1mean]);
% % % scatter([1*ones(1,length(amygrows)),2*ones(1,length(lpporrows)),3*ones(1,length(dlgnporrows)),4*ones(1,length(lpv1rows)),5*ones(1,length(dlgnv1rows))],...
% % %     [amygvals,lpporvals,dlgnporvals,lpv1vals,dlgnv1vals]);
% % % ylim([0.1 4.5])
% % % hold all
% % % plot([0.5 npaths+0.5],[1 1],'r--')

close all
bar([amygmean,lppormean,dlgnpormean,lplimean,dlgnlimean,lplmmean,dlgnlmmean,lpv1mean,dlgnv1mean],'BaseValue',1)

ylabel('Patch/interpatch OD ratio')
set(gca,'XTick',1:npaths)
xlim([0.5 npaths+0.5])
ylim([0.3 4.5])
% set(gca,'XTickLabel',{'LA-->POR','LP-->POR','dLGN-->POR','LP-->V1','dLGN-->V1','V1-->POR'})
set(gca,'XTickLabel',{'LA-->POR','LP-->POR','dLGN-->POR','LP-->LI','dLGN-->LI','LP-->LM','dLGN-->LM','LP-->V1','dLGN-->V1'});
 set(gca,'YScale','log')
 set(gca,'YTick',[0.5 1 2 4])
 hold all
 errorbar(1:npaths,[amygmean,lppormean,dlgnpormean,lplimean,dlgnlimean,lplmmean,dlgnlmmean,lpv1mean,dlgnv1mean],...
    [amygster,lpporster,dlgnporster,lplister,dlgnlister,lplmster,dlgnlmster,lpv1ster,dlgnv1ster],'r','LineStyle','none')