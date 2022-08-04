%-----------------------------------------------------------------------------------------------------------------------
%-- CTResponseCurve.m -- Plots a contrast response curve,for sine
%and square waves
%--	GCD, 1/12/06
%-----------------------------------------------------------------------------------------------------------------------
function CTResponseCurve(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;		%needed for defines like IN_T1_WIN_CD
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01
Path_Defs;

symbols = {'ko' 'k*' 'go' 'mo' 'b*' 'r*' 'g*' 'c*'};
lines = {'k-' 'k--' 'g-' 'm-' 'b--' 'r--' 'g--' 'c--'};


%get the column of values of contrasts in gratings_params matrix
contrast = data.gratings_params(GRAT_CONTRAST,:,PATCH1);

%get indices of any NULL conditions (for measuring spontaneous activity)
null_trials = logical( (contrast == data.one_time_params(NULL_VALUE)) );

%get the column of stimulus type values
stim_type = data.gratings_params(GRAT_TYPE,:,PATCH1);
unique_stim_type = munique(stim_type(~null_trials)');

%now, get the firing rates for all the trials 
spike_rates = data.spike_rates(SpikeChan, :);

unique_contrast = munique(contrast(~null_trials)')

%now, remove trials that do not fall between BegTrial and EndTrial
trials = 1:length(contrast);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

% Calculate spontaneous rates before looping through so can calculate DTI
null_resp = data.spike_rates(SpikeChan, null_trials & select_trials);
null_rate = mean(null_resp);

figure(2);
set(gcf,'PaperPosition', [.51 .51 20 27.5], 'Position', [450 50 500 573], 'Name', 'Grating Contrast Response Curve');
subplot(2, 1, 2);

for j = 1:length(unique_stim_type)

    grating_select = logical( (stim_type == unique_stim_type(j)) );
    
    plot_x=contrast(~null_trials & select_trials & grating_select);
    plot_y=spike_rates(~null_trials & select_trials & grating_select);

    %NOTE:  Inputs to PlotTuningCurve.m must be column vectors, not row vectors
    %because of munique().
    [px, py, perr] = PlotTuningCurve(plot_x', plot_y', symbols{j}, lines{j}, 1, 0);
    %NOTE: last argument=0 means just get output, no plot

    unique_contrasts = px; %store category groups for ANOVAs

    %Plot the data
    hold on;
    plot(plot_x, plot_y, 'k.');
    errorbar(px, py, perr, perr, symbols{j});
    hold on;

    %compute discrimination index
    [CDI(j), var_term(j)] = Compute_DDI(plot_x, plot_y);
    
    %Now, fit the data with a hyperbolic ratio function and plot this as
    %well
    means = [px py];
    raw = [plot_x' plot_y'];
    [pars] = hyperbolicratio_fit(means, raw);
    x_interp = [(px(1)):0.01:(px(length(px)))];
    y_interp = hyperbolicratio_func(x_interp, pars);
    plot(x_interp, y_interp, 'k-');
    
   
    %now, get the firing rate for NULL condition trials and add spontaneous
    %rate to plot
    null_x = [min(px) max(px)];
    null_resp = data.spike_rates(SpikeChan, null_trials & select_trials);
    null_rate = mean(null_resp);
    null_y = [null_rate null_rate];
    hold on;
    plot(null_x, null_y, 'k--');
    hold off;

    yl=ylim;
    ylim([0 yl(2)]); %set lower limit of the Y axis to zero
    xlabel('Grating Contrast');
    ylabel('Response (spikes/sec)');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Added on 04/21/06, suppressed until SF_loggaussfit is working
    %Compute R^2 for both means and raw values
    y_fit = hyperbolicratio_func(px, pars);
    y_fit(y_fit<0) = 0;
    %add a column of ones to yfit to make regress happy
    y_fit = [ones(length(y_fit), 1) y_fit];
    [b, bint, r, rint, stats1(j,:)] = regress(py, y_fit);

    y_fit_raw = hyperbolicratio_func(plot_x', pars);
    y_fit_raw(y_fit_raw<0) = 0;
    y_fit_raw = [ones(length(y_fit_raw), 1) y_fit_raw];
    [b, bint, r, rint, stats2(j,:)] = regress(plot_y', y_fit_raw);

    %Calculate some metrics and stats to print in the plot

    p_value(j) = spk_anova(plot_y, plot_x, unique_contrasts);
    Rmax(j) = pars(1);
    C50(j) = pars(2);
    n(j) = pars(3);
    DC_offset(j) = pars(4);
    max_rate(j) = max(y_interp);

    CMI(j) = 1 - (DC_offset - null_rate)/(max_rate - null_rate);

    %Do chi-square goodness of fit test
    [chi2(j), chiP(j)] = Chi2_Test(plot_x, plot_y, 'hyperbolicratio_func', pars, length(pars));

    [pmax{j}.y max_i] = max(y_interp);
    pmax{j}.x = x_interp(max_i);
    [pmin{j}.y min_i] = min(y_interp);
    pmin{j}.x = x_interp(min_i);
    
    [pmax{j}.y max_i] = max(y_interp);
    pmax{j}.x = x_interp(max_i);
    [pmin{j}.y min_i] = min(y_interp);
    pmin{j}.x = x_interp(min_i);

    index = 1;
    ysearch = y_interp(index);
    while (ysearch < (pmin{j}.y + (pmax{j}.y-pmin{j}.y)/2))
        ysearch = y_interp(index);
        if (index == length(y_interp))
            break;
        else
            index = index + 1;
        end
    end
    halfmax_slope(j) = (ysearch - pmin{j}.y)/(x_interp(index)-pmin{j}.x);

end

%now, print out some useful information in the upper subplot
subplot(2, 1, 1);
PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

%now print out useful values
% pmax, pmin, py,
PrintGratingContrastData(p_value, Rmax, C50, n, DC_offset, max_rate, stats1, stats2, chi2, chiP, CDI, CMI);

output = 1;
if (output == 1)

    %------------------------------------------------------------------------
    %write out all relevant parameters to a cumulative text file, GCD 8/08/01
    %write out one line for each stimu_type for each neuron.
    outfile = [BASE_PATH 'ProtocolSpecific\GratingContrast\GratingContrast_Summary.dat'];
    printflag = 0;
    if (exist(outfile, 'file') == 0)    %file does not yet exist
        printflag = 1;
    end
    fid = fopen(outfile, 'a');
    if (printflag)
        fprintf(fid, 'FILE\t\t\t PrfOR\t PrfSF\t PrfTF\t RFXctr\t RFYctr\t RFDiam\t GrType\t MaxLoc\t MaxRsp\t MinLoc\t MinRsp\t Spont\t CMI\t\t AnovaP\t SFDI  \t VarTerm\t MaxX \t MinX \t RsqMn\t PMean\t RsqRaw\t Praw\t\t Chi2\t\t ChiP\t\t Rmax\t\t\t C50\t\t n\t\t Offset\t Halfmax_slope\t');
        fprintf(fid, '\r\n');
        printflag = 0;
    end
    for j = 1:length(unique_stim_type)
        buff = sprintf('%s\t %5.1f\t %5.2f\t %5.3f\t %5.2f\t %5.2f\t %5.2f\t %5.1f\t %5.4f\t %5.2f\t %5.4f\t %5.2f\t %5.2f\t %5.3f\t %10.8f\t %5.3f\t %6.3f\t %5.2f\t %5.2f\t %7.4f\t %6.4f\t %10.8f\t %6.4f\t %10.8f\t %6.4f\t %10.8f\t %5.2f\t %5.2f\t %6.4f\t %6.4f\t %6.4f\t %6.4\t', ...
            FILE, data.neuron_params(PREFERRED_ORIENTATION, 1), data.neuron_params(PREFERRED_SPATIAL_FREQ, 1), data.neuron_params(PREFERRED_TEMPORAL_FREQ, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
            unique_stim_type(j), pmax{j}.x, pmax{j}.y, pmin{j}.x, pmin{j}.y, null_rate, CMI(j), p_value(j), CDI(j), var_term(j), max(unique_contrast), min(unique_contrast), stats1(j,1), stats1(j,3), stats2(j,1), stats2(j,3), chi2(j), chiP(j), Rmax(j), C50(j), n(j), DC_offset(j), halfmax_slope(j) );
        fprintf(fid, '%s', buff);
        fprintf(fid, '\r\n');
    end
    fclose(fid);
    %------------------------------------------------------------------------

end  %if (output == 1)

return;


