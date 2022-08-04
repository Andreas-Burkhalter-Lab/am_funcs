%%% get corr results for multiple image pairs

function res = tabulate_cor(secnums,images1,images2,rois,filtsizes)

rois_is_file_basename = 1;

nimages = length(images1);
nfiltsizes = length(filtsizes);

if rois_is_file_basename
    for i = 1:nimages
        roinames{i} = [rois{1} num2str(secnums(i)) rois{2}];
    end
else
    roinames = rois;
end

if size(secnums,2)>1
    secnums = secnums';
end

filtvarnames = cellfun(@(x)['cor_filt_' x], cellstr(num2str([filtsizes]'))','UniformOutput',false);
res = table(secnums,'VariableNames',{'sec'});
res{:,2:nfiltsizes+2} = NaN(nimages,nfiltsizes+1);
res.Properties.VariableNames(2) = {'cor_unfilt'};
res.Properties.VariableNames(3:end) = strrep(filtvarnames,' ','');

for i = 1:nimages
    c = pixel_correlation(images1{i},images2{i},roinames{i},filtsizes);
    res{i,2:nfiltsizes+2} = c.corrTable.corrCoef';
end