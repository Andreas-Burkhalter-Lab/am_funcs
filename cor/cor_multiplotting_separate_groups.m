%%% cor_multiplotting_separate_groups: plot 4 axes without overlaying p-p with ip-ip
% may no longer be functional
%
%%%% updated 2020/8/28

corfig = figure('units','normalized','outerposition',pars.fig_outerposition,'WindowStyle',windowstyle); 
suptitle(titlestring)
clear cortable_out

subplot(subplot_rowcol(1), subplot_rowcol(2), 1)
pars.var_to_plot = 'corcoef_rest'; % dff cor during dark period
pars.seedquant = [1 2 3]; % interpatches
[cortable_out{1},data_out{1}] = pairwise_cor_plotting(cortable,pars);
ylim(pars.rest_ylim)
title('interpatch rest')

subplot(subplot_rowcol(1), subplot_rowcol(2), 2)
pars.var_to_plot = 'corcoef_noise'; % noise cor during stim
pars.seedquant = [1 2 3]; % interpatches
pairwise_cor_plotting(cortable,pars);
ylim(pars.noise_ylim)
title('interpatch noise')

subplot(subplot_rowcol(1), subplot_rowcol(2), 3)
pars.var_to_plot = 'corcoef_rest'; % dff cor during dark period
pars.seedquant = [4 5 6]; % patches
[cortable_out{2},data_out{2}] = pairwise_cor_plotting(cortable,pars);
ylim(pars.rest_ylim)
title('patch rest')

subplot(subplot_rowcol(1), subplot_rowcol(2), 4)
pars.var_to_plot = 'corcoef_noise'; % noise cor during stim
pars.seedquant = [4 5 6]; % patches
pairwise_cor_plotting(cortable,pars);
ylim(pars.noise_ylim)
title('patch noise')