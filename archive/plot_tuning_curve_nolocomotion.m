%%% plot tuning curve of selected parameter of  selected ROI from tuning table
%
% plot_tuning_curve(tuning_table,tablerow,stimparname)
% stimparname = 'sf','tf', or 'orient'
% updated 18/05/11 on thermaltake


function plot_tuning_curve(tuning_table,tablerow,stimparname)

nxvals = 1e5;
new_fig = 1;

trow = tuning_table(tablerow,:);
tdat = trow{1,[stimparname '_trials']}{1};
stimparvals = tdat{:,stimparname};
xvals = linspace(min(stimparvals),max(stimparvals),nxvals);
resp = tdat{:,'resp'};
ntrials = size(resp,2);

if new_fig
    figure
end

switch stimparname
    case 'sf'
        fitfunc = @(q,stimval)q(2).*exp(-1./2./q(4).^2.*(log((stimval+q(5))./(q(3)+q(5)))).^2)+q(1); % lognormal from Gao et al. 2010
        xlabel('Spatial frequency (cyc/deg)')
        xscaletype = 'log';     
    case 'tf'
        fitfunc = @(q,stimval)q(2).*exp(-1./2./q(4).^2.*(log((stimval+q(5))./(q(3)+q(5)))).^2)+q(1); % lognormal from Gao et al. 2010
        xlabel('Temporal frequency (hz)')
        xscaletype = 'log';     
    case 'orient'
        fitfunc = @(q,stimval)q(1).* exp( q(2).* (cos(stimval.*pi/180-q(3)*pi/180)-1)) + q(4).* exp( q(5).* (cos(stimval.*pi/180-q(6)*pi/180)-1)) + q(7);  % von mises from Gao et al. 2010                    
        xlabel('Orientation (deg)')
        xscaletype = 'linear';     
    otherwise 
        error('unknown stim par name')
end

plot(xvals,fitfunc(trow{1,[stimparname '_fitparams']}{1},xvals),'r')
hold on
title( ['row ' num2str(tablerow) ', anovap = ' num2str(trow{1,[stimparname '_anovap']})] )
set(gca, 'XScale', xscaletype)
datplot = [repmat(stimparvals,ntrials,1), resp(:)];
datplot( isnan(datplot(:,2)),: ) = []; % delete nan trials

scatter(datplot(:,1),datplot(:,2))
ylabel('Response dF/F')
hold off