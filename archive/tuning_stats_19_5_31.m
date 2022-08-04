%% this version does parameter definition in a messier way than future versions

%%%%%%%%%%% perform anova and correlation tests on gcamp rois
%%%% for each table row, at least one criteria_var from each criterai var list must be true for a table row to be analyzed
%%%% varname == the test statistic
%%%% sorting_var == variable used to split rois into groups
%%%%%%%%%% updated on thermaltake 2019-05-03

% close all
clear meets_criteria varname criteria_var1 criteria_var2 criteria_var3 criteria_var_quantile subject % other_criteria

plot_bar_with_sem = 0; 
plot_anova = 'off';

% varname = 'tf_pref';
% varname = 'tf_hwhm';
% varname = 'tf_anovap';
% criteria_var1 = 'tf_sgnf';


% varname = 'sf_pref';
% varname = 'sf_hwhm';
% varname = 'sf_anovap';
% criteria_var1 = 'sf_sgnf';

% varname = 'orient_hwhm';
% varname = 'orient_sgnf';
% varname = 'orient_anovap';
% criteria_var1 = 'orient_sgnf';

% varname = 'locm_cor_r';
varname = 'pakan_loc_index';
% varname = 'stim_overlap_frac';
% varname = 'stim_on_rf';
% varname = 'pupilwidth_cor_r'; 
% varname = 'pupilwidth_cor_p';
% criteria_var1 = 'pupilwidth_cor_sgnf';
 
criteria_var2 = {'sf_sgnf','tf_sgnf','orient_sgnf','rspv'};
% criteria_var2 = {'sf_sgnf','tf_sgnf','orient_sgnf','rf_sgnf'}; %%% if stimulus is tuned for any stimulus parameters we tested
% criteria_var2 = 'rspv';
% criteria_var2 = 'stim_on_rf';
criteria_var3 = 'rf_sgnf'; 
% varname = 'rspv';

% other_criteria = roitable.locm_cor_p > 0.05; 
% other_criteria = roitable.stim_overlap_frac > 0.0; 
% other_criteria = roitable.rf_overlap_frac_aprx > 0.1;
% other_criteria = roitable.stim_overlap_frac == 0; 


plane = [4 3  2 1]; 
layer = [23];
% other_criteria = strcmp(roitable.day, filetable.day{1});
% subject = [18164]; 
% subject = [18165]; 
% subject = [19002]; 
% subject = [19042]; 
criteria_var_quantile = 'quantile'; %%% don't count cells that fall in areas excluded from patchiness analysis, which are coded as quantile==0
clear meets_criteria
criteria_var1 = vardefault('criteria_var1','sub'); % if criteria var isn't selected, don't filter out any rois
criteria_var2 = vardefault('criteria_var2','sub'); % if criteria var isn't selected, don't filter out any rois
criteria_var3 = vardefault('criteria_var3','sub'); % if criteria var isn't selected, don't filter out any rois
other_criteria = vardefault('other_criteria',true(height(roitable),1));
subject = vardefault('subject',unique(roitable.sub))'; %%% if no subject specified, analyze all subjects
meets_criteria = any(roitable{:,criteria_var1} , 2) & ...
                 any(roitable{:,criteria_var2} , 2) & ...
                 any(roitable{:,criteria_var3} , 2) & ...
                 other_criteria & ...
             logical(roitable{:,criteria_var_quantile}) ; 

meets_criteria = meets_criteria & any(logical(roitable.plane == plane) , 2); %% match any of the specified planes
meets_criteria = meets_criteria & any(logical(roitable.layer == layer) , 2); %% match any of the specified layers
meets_criteria = meets_criteria & any(logical(roitable.sub == subject) , 2); %% match any of the specified layers
meets_criteria = meets_criteria & ~isnan(roitable{:,varname}); %% exclude rois that have value NaN for the variable of interest
n_rois_meet_criteria = nnz(meets_criteria);
% sorting_var = 'plane';
sorting_var = 'inpatch';
% sorting_var = 'quantile';

% get values for each group so standard error can be calculated
labelnames = unique(roitable{meets_criteria,sorting_var});
sorted_vals = table(labelnames, cell(size(labelnames)), NaN(size(labelnames)), NaN(size(labelnames)), 'VariableNames', {sorting_var, varname, 'mean', 'sem'});
for ilabel = 1:length(labelnames)
    thesevals = roitable{meets_criteria & roitable{:,sorting_var}==labelnames(ilabel),varname};
    sorted_vals{ilabel,varname} = {thesevals};
    sorted_vals.mean(ilabel) = mean(thesevals);
    sorted_vals.sem(ilabel) = std(thesevals) ./ sqrt(length(thesevals));
end
% do anova
values_to_test = double(roitable{meets_criteria,varname}); 
labels_to_test = double(roitable{meets_criteria,sorting_var});
[anovap, anovatab, anovastats] = anova1(values_to_test, labels_to_test, plot_anova)

[r p] = corrcoef(values_to_test, labels_to_test)

if plot_bar_with_sem
    figure
    bar(sorted_vals.mean)
    hold on
    errorbar(1:height(sorted_vals),sorted_vals.mean,sorted_vals.sem,'.')
    hold off
end

