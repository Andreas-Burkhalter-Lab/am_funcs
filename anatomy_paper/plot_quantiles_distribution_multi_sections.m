%%% plot cell density in quantiles for each section listed in the filetable

function plot_quantiles_distribution_multi_sections(filetable)

rows = 3;
columns = 4;

nsections = height(filetable);
figure
for i = 1:nsections
    subplot(rows,columns,i)
    bar(1:height(filetable.quantcells{i}),filetable.quantcells{i}.cells_per_sqmm)
    title(['sec', num2str(filetable.sec(i)), 'chi_p=', num2str(filetable.chi_p(i))])
    ylabel('cells_per_mm^2')
end