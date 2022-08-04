%%% plot tuning curve of selected parameter of  selected ROI from tuning table
%
% plot_tuning_curve(tuningdata,tablerow,stimparname)
% stimparname = 'sf','tf', or 'orient'
% updated 20/9/23 on thermaltake


function tuning_results = plot_tuning_curve(tuningdata, tablerow, stimparname, locomotion_forw_thresh_mps, pars)

if ~exist('pars','var') %%% using vardefault.m here sometimes crashes Matlab; related to power usage?
    pars = struct;
end
pars.nxvals = field_default(pars,'nxvals',1e5);
pars.new_fig = field_default(pars,'show_plot',1); % if false, do not plot any data, only output the fit results
pars.new_fig = field_default(pars,'new_fig',0);
pars.get_hwhm = field_default(pars,'get_hwhm',0); %% if true, compute half width at half maximum to output in tuning_results
pars.show_roi_location = field_default(pars,'show_roi_location',0); 
pars.LineWidth = field_default(pars,'LineWidth',2); % fitted curve line width
pars.Marker = field_default(pars,'Marker','.'); % marker type of trial points
pars.SizeData = field_default(pars,'SizeData',200); % size of trial points
% if not specified from input arg, set locomotion threshold defaults for dividing high- vs. low-locomotion trials; value from Pakan et al. 2016 is 0.001
    default_locomotion_forw_thresh_mps_ablab = 0.001; % mean forward locomotion over the course of the trial in meters per second required to plot a trial as high locomotion
    default_locomotion_forw_thresh_mps_hanlab = 0.008; 
pars.low_locm_plotcolor = field_default(pars,'low_locm_plotcolor',[1 0 0]);
pars.high_locm_plotcolor = field_default(pars,'high_locm_plotcolor',[0 0 0]);
pars.show_title = field_default(pars,'show_title',1);
    
try tuningdata.stimpars = tuningdata.tuningpars.stimpars; end
try tuningdata.locm_trials = tuningdata.tuningpars.locm_trials; end
if ~exist('locomotion_forw_thresh_mps','var') || isempty(locomotion_forw_thresh_mps) % if locomotion threshold not specified, use the default
    if strcmp(tuningdata.stimpars.computer,'IMAGING_VR-PC')
        locomotion_forw_thresh_mps = default_locomotion_forw_thresh_mps_hanlab;
    elseif strcmp(tuningdata.stimpars.computer,'ANDREWLAB-PC')
        locomotion_forw_thresh_mps = default_locomotion_forw_thresh_mps_ablab;
    end
end

trow = tuningdata.tuning(tablerow,:);
tdat = trow{1,[stimparname '_trials']}{1};
stimparvals = tdat{:,stimparname};
pars.xvals = linspace(min(stimparvals),max(stimparvals),pars.nxvals);
resp = tdat{:,'resp'};
ntrials = size(resp,2);
locm_trials = tuningdata.locm_trials.(stimparname);
locm_superthresh = locm_trials.locm_forw_mps > locomotion_forw_thresh_mps; 
stimparval_matrix = repmat(locm_trials{:,stimparname},1,tuningdata.stimpars.repetitions);

% low locomotion
resp_low_locm = resp; 
resp_low_locm(locm_superthresh) = NaN; % include only subthresh trials
meanforfit_low_locm(:,1) = stimparvals; % first col = stimparval, second = response val
meanforfit_low_locm(:,2) = nanmean(resp_low_locm,2); % first col = stimparval, second = response val
respforfit_low_locm = NaN(nnz(~locm_superthresh),2); % stimpars and responses for low locomotion trials
respforfit_low_locm(:,1) = stimparval_matrix(~locm_superthresh);  % first col = stimparval, second = response val
respforfit_low_locm(:,2) = tdat.resp(~locm_superthresh); % first col = stimparval, second = response val
fitresults_low_locm = get_function_fit(respforfit_low_locm, meanforfit_low_locm, stimparname, pars.get_hwhm);
if pars.show_plot
    if pars.new_fig
        figure
    end
    plothandle_low_locm = make_fit_plot(fitresults_low_locm, tablerow, pars, pars.high_locm_plotcolor); % plot low loc trials    
end
% high locomotion
resp_high_locm = resp; 
resp_high_locm(~locm_superthresh) = NaN; % include only subthresh trials
meanforfit_high_locm(:,1) = stimparvals; % first col = stimparval, second = response val
meanforfit_high_locm(:,2) = nanmean(resp_high_locm,2); % first col = stimparval, second = response val
if nnz(~isnan(meanforfit_high_locm(:,2))) > 2 % if we have high locomotion trials for at least 2 different stim param values, function fit the high locomotion trials
    respforfit_high_locm = NaN(nnz(locm_superthresh),2);
    respforfit_high_locm(:,1) = stimparval_matrix(locm_superthresh);  % first col = stimparval, second = response val
    respforfit_high_locm(:,2) = tdat.resp(locm_superthresh); % first col = stimparval, second = response val
    fitresults_high_locm = get_function_fit(respforfit_high_locm, meanforfit_high_locm, stimparname, pars.get_hwhm);
    if pars.show_plot
        hold on
        plothandle_high_locm = make_fit_plot(fitresults_high_locm, tablerow, pars, pars.low_locm_plotcolor); % plot high loc trials
    %     legend([plothandle_low_locm, plothandle_high_locm], {'low locomotion','high locomotion',},'Position',[0.23 0.7987 0.1184 0.1005])
        legend([plothandle_low_locm, plothandle_high_locm], {'low locomotion','high locomotion',})
    end
    tuning_results.fitresults_high_locm = fitresults_high_locm; 
else
    tuning_results.fitresults_high_locm = struct; 
end

%%% output results 
tuning_results.fitresults_low_locm = fitresults_low_locm; 


%%% fit tuning curves to the trial response data
function fitout = get_function_fit(respforfit, meanforfit, stimparname, get_hwhm)
    respforfit( isnan(respforfit(:,2)),: ) = []; % delete nan trials before fitting
    fitout.stimparname = stimparname;
    fitout.respforfit = respforfit;
    switch stimparname
        case 'sf'
            fitparams = SF_loggaussfit(meanforfit,respforfit);
            fitout.sf_pref = fitparams(3); 
            fitout.sf_fit_base_F = fitparams(1);
            fitout.sf_fit_amp = fitparams(2); % amplitude
            fitout.sf_fit_sigma = fitparams(4); % standard deviation
            fitout.sf_fit_log_offset = fitparams(5);
            fitout.prefval = fitout.sf_pref;
            fitout.fitfunc = @(q,stimval)q(2).*exp(-1./2./q(4).^2.*(log((stimval+q(5))./(q(3)+q(5)))).^2)+q(1); % lognormal from Gao et al. 2010
            fitout.hwhm_xspacing = 0.0005;
            fitout.hwhm_xmax = 3.5; 
            fitout.xlabeltext = ('Spatial frequency (cyc/deg)'); 
%             fitout.xscaletype = 'linear'; 
            fitout.xscaletype = 'log'; 
        case 'tf'
            fitparams = TF_loggaussfit(meanforfit,respforfit);
            fitout.tf_pref = fitparams(3);
            fitout.tf_fit_base_F = fitparams(1);
            fitout.tf_fit_amp = fitparams(2); % amplitude
            fitout.tf_fit_sigma = fitparams(4); % standard deviation
            fitout.tf_fit_log_offset = fitparams(5);
            fitout.prefval = fitout.tf_pref;
            fitout.fitfunc = @(q,stimval)q(2).*exp(-1./2./q(4).^2.*(log((stimval+q(5))./(q(3)+q(5)))).^2)+q(1); % lognormal from Gao et al. 2010
            fitout.hwhm_xspacing = 0.01; 
            fitout.hwhm_xmax = 30;  
            fitout.xlabeltext = ('Temporal frequency (hz)');
            fitout.xscaletype = 'log'; 
        case 'orient'
            fitparams = von_mises_fit(meanforfit,respforfit);
            fitout.fitfunc = @(q,stimval)q(1).* exp( q(2).* (cos(stimval.*pi/180-q(3)*pi/180)-1)) + q(4).* exp( q(5).* (cos(stimval.*pi/180-q(6)*pi/180)-1)) + q(7);  % von mises from Gao et al. 2010                    
            if fitout.fitfunc(fitparams,fitparams(3)) < fitout.fitfunc(fitparams,fitparams(6)) % if the second pref orient generates higher response than first pref orient, switch their order
                fitparams = [fitparams(4:6); fitparams(1:3); fitparams(7)]; % put most preferred orient first
            end
            fitout.orient_pref1 = fitparams(3);
            fitout.orient_amp1 = fitparams(1);
            fitout.orient_width1 = fitparams(2);
            fitout.orient_pref2 = fitparams(6);
            fitout.orient_amp2 = fitparams(4);
            fitout.orient_width2 = fitparams(5);
            fitout.orient_fit_base_F = fitparams(7);
            fitout.prefval = fitout.orient_pref1;
            fitout.hwhm_xspacing = 0.1;
            fitout.hwhm_xmax = fitout.orient_pref1 + 200;
            fitout.xlabeltext = ('Orientation (deg)');
            fitout.xscaletype = 'linear'; 
        case 'diam'
            fitparams = diff_erf_fit(meanforfit,respforfit);  %%% fitting function from Gao et al. 2011
            fitout.fitfunc = @(q,stimval)q(1)*erf(stimval/q(2))-q(3)*erf(stimval/(q(2)+q(4)))+q(5);     %%% difference-of-error function plus baseline
            fitout.diam_excit_size = fitparams(2); % size of the excitatory center
            fitout.diam_excit_amp = fitparams(1);
            fitout.diam_inhib_size = fitparams(4);% size of the inhibitory surround
            fitout.diam_inhib_amp = fitparams(3);
            fitout.diam_fit_base_F = fitparams(5);
            fitout.Properties.UserData = 'stimulus diameter in degrees';
            fitout.prefval = fitout.diam_excit_size;
            fitout.hwhm_xspacing = 0.05; 
            fitout.hwhm_xmax = 100; 
            fitout.xlabeltext = ('Diameter (deg)');
            fitout.xscaletype = 'linear';
        otherwise 
            error('unknown stim par name')
    end
    if get_hwhm
        if fitout.prefval+fitout.hwhm_xspacing >= fitout.hwhm_xmax
            fitout.hwhm = NaN; % tuning curve is so broad, don't bother computing hwhm
        else
            fitout.hwhm = hwhm(fitout.fitfunc, fitparams, [fitout.prefval  fitout.hwhm_xmax], fitout.hwhm_xspacing, fitout.prefval, fitout.([stimparname,'_fit_base_F'])); 
        end
    end
    
    fitout.fitparams = fitparams;
end

%%% show location of roi within time-avg image
if pars.show_roi_location 
    figure
    if strcmp(tuningdata.stimpars.computer,'IMAGING_VR-PC') % de-rotate
        ee = edge(full(trow.roi_image_pre_rotate{:}));
    elseif strcmp(tuningdata.stimpars.computer,'ANDREWLAB-PC') % no de-rotation
        ee = edge(full(trow.roi_image_prereg{:}));
    end
    [y x] = find(ee);
    imagesc(tuningdata.meanImage_pre_rotate)
    colormap gray
    hold on
    scatter(x,y,'r','.')
    hold off
end



%%%% plot data points and tuning curve
function [plothandle] = make_fit_plot(fitresults, tablerow, pars, plotcolor)
    plothandle = plot(pars.xvals,fitresults.fitfunc(fitresults.fitparams,pars.xvals),'Color',plotcolor);
    plothandle.LineWidth = pars.LineWidth;
    hold on
    if pars.show_title
        title( [fitresults.stimparname,' , row ' num2str(tablerow) ', anovap = ' num2str(trow{1,[fitresults.stimparname '_anovap']})] ) % title with anova p
    %     title( ['row ', num2str(tablerow), ', ', fitresults.stimparname, '_tuning'])
    end
    set(gca, 'XScale', fitresults.xscaletype)
    scplot = scatter(fitresults.respforfit(:,1), fitresults.respforfit(:,2), 'CData', plotcolor);
    scplot.Marker = pars.Marker; 
    scplot.SizeData = pars.SizeData;
    scplot.MarkerFaceColor = plotcolor;
    ylabel('Response dF/F')
    xlabel(fitresults.xlabeltext)
    hold off
end

end