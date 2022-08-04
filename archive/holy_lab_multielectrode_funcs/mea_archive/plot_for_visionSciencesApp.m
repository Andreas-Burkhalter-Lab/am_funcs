%%% plot sizetuning for vision sciences application submitted mid-October
%%% 2015

% set(0,'DefaultFigureWindowStyle','docked');
% cd /mnt/andrew_vivid/recordings/08152015/1sizetuning
% m = merecmm('1sizetuning.merec');
% load oanalyzed_ssstimdata15-Aug-2015_1sizetuning.mat
% close all
% hold off
clf

%%% divide response spikes by 2 to get in hz (stim dur and spike summation
%%%     time = 2s)
%%% uX = unit index in unitAvg_tffixed to plot
u1 = 15;% 4 is not surround-suppressed
u2 = 4; % 15 is surround suppressed
u1hz = 0.5.*double(unitAvg_tffixed{u1,1});
u2hz = 0.5.*double(unitAvg_tffixed{u2,1});

% gaussian fitting
u1gauss1 = fit(diam_vals',u1hz','gauss1');
u2gauss1 = fit(diam_vals',u2hz','gauss1');
u1gauss2 = fit(diam_vals',u1hz','gauss2');
u2gauss2 = fit(diam_vals',u2hz','gauss2');

% get spikes per trial for errorbars
u1trialsps = NaN(pars.repetitions,pars.n_diams);
u2trialsps = NaN(pars.repetitions,pars.n_diams);
for diamind = 1:pars.n_diams
    u1trialsps(:,diamind) = 0.5 .* [angleAvg_tffixed.sps_all{diamind}(:,u1)]'; % spikes on each trial
    u2trialsps(:,diamind) = 0.5 .* [angleAvg_tffixed.sps_all{diamind}(:,u2)]'; % spikes on each trial
end
u1sem = std(u1trialsps) ./ sqrt(pars.n_diams);
u2sem = std(u2trialsps) ./ sqrt(pars.n_diams);
u1means = mean(u1trialsps); % should be the same as u1hz
u2means = mean(u2trialsps); % should be the same as u2hz

% plotting
hold all
markersize = 10;
ebarswidth = 0.1; % default width (0.35) seems to be minimum
gausscurves_width = 1.7; 
axis_font_size = 35;
tickfontsize = 23;
legfontsize = 11;
errorbar(diam_vals,u1means,u1sem,u1sem,'.k','MarkerSize',markersize,'LineWidth',ebarswidth)
errorbar(diam_vals,u2means,u2sem,u2sem,'.r','MarkerSize',markersize,'LineWidth',ebarswidth)
u1gauss_ax = plot(u1gauss1,'k');
u2gauss_ax = plot(u2gauss1,'r');
set(u1gauss_ax,'LineWidth',gausscurves_width)
set(u2gauss_ax,'LineWidth',gausscurves_width)
% plot(diam_vals,u1hz) % plot diam on linear scale
% plot(diam_vals,u2hz) % plot diam on linear scale
% plot(unitAvg_tffixed{u1,1}) % plot diam on log scale (original distribution of diam vals) 
% plot(unitAvg_tffixed{u2,1}) % plot diam on log scale (original distribution of diam vals)
xlabel('Stimulus diameter (deg)','FontSize',axis_font_size)
ylabel('Response (spikes/s)','FontSize',axis_font_size)
set(gca,'FontSize',tickfontsize)

xlim([2 52])

% legend({'Strong surround-suppresion','Weak surround suppresion'})
% leg = legend([u1gauss_ax u2gauss_ax],...
%     {'Surround-suppressed','Not surround-suppressed'},'FontSize',legfontsize);
legend off

% hold all
% for i = 1:16
%     plot(unitAvg_tffixed{i,1})
% end
