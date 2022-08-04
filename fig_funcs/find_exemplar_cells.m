%%%% for finding examplar cells (for locomotion modulation fig)
% first load cor20-8-6.mat, run tuning stats, require sf/tf/orient_sgnf and quant == 1-2, 3-4, or 5-6
%   or require sf/tf/orient p < 0.01 for better tuning curves
%
%%% input these cells to exemplar_cells.xlsx then run locmod_tuning_curves to plot the cells
%
% updated 2020/8/18

%%%%% variable to search for
tuning_var = 'sf';
% tuning_var = 'tf';
% tuning_var = 'orient';

rr = roitable(ops.meets_criteria,:); % keep only potential examplar cells
[~, ordr] = sort(rr.pakan_loc_index); % sort by locomotion modulation index
rr = rr(flipud(ordr),:);
rr.day = datetime(rr.day,'Format','yyyy-MM-dd'); % make easier to copy into excel
rr = movevars(rr,'quantile','Before','sub');
rr = movevars(rr,'centeryx_prereg','After','plane');
rr = movevars(rr,{[tuning_var, '_anovap'], 'pakan_loc_index'},'After','centeryx_prereg');