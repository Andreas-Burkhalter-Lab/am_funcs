%%%% last updated 5/11/18

cell_index = 5;
show_anova_plots = 'on';
for cell_index = 12%1:height(p3.tuning)
cell_index
    [combo_table IA combo_index] = unique(p3.stimpar_sets(:,{'Angle','diam','sfreq','tfreq'})); % unique param combinations
    ncombos = height(combo_table);
    combo_table.comboname =  [1:ncombos]';
    combo_table.resp_mean = NaN(ncombos,1);
    combo_table.resp = NaN(ncombos,p3.stimpars.repetitions);
    combonames_for_anova = kron(combo_table.comboname,ones(size(combo_table.resp,2),1));

    for icombo = 1:ncombos
        combo_table.resp(icombo,:) = [p3.stimpar_sets.F_during_stim(find(combo_index==icombo),cell_index)]';
        combo_table.resp_mean(icombo,1) = mean(combo_table.resp(icombo,:)); 
    end

    resp_for_anova = combo_table.resp'; % transpose so that resps will be extracted in order of combo type, not trial number
    resp_for_anova = resp_for_anova(:);
    %%% add baseline data - labeled as zeros
    resp_for_anova = [resp_for_anova; p3.stimpar_sets.F_prestim(:,cell_index)]; % append baseline
    combonames_for_anova = [combonames_for_anova; zeros(height(p3.stimpar_sets),1)]; % append baseline labels

    close all
    [p,tbl,stats] = anova1(resp_for_anova, combonames_for_anova, show_anova_plots);
    p3.tuning.rspP(cell_index) = p;
    p3.tuning.rsp1(cell_index) = p<.05;
end
