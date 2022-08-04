%-----------------------------------------------------------------------------------------------------------------------
%-- SizeTuningCurve.m -- Plots a size tuning curve,for sine
%and square waves
%--	CMA, 5/12/06
%-----------------------------------------------------------------------------------------------------------------------

%%% note: this script automatically looks for 'grating sizes,' not 'dot
%%%    size'; for dot size tuning, use DotSizeTuningCurve.m instead


%%%% 5/19/16 AM added conditional when writing outfile to write into the correct
%%%%    directory if on a specific computer; functionality should be
%%%% 5/19/16 AM edited formatting of output file columns
%%% 3/14/16 AM corrected capitalization errors
%%%% 8/23/16 AM added optional 'skip_trials' argin - these trials will not be analyzed;
%%%%    trials not included within BegTrial-EndTrial will still be skipped

function SizeTuningCurve(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, skip_trials)

TEMPO_Defs;		%needed for defines like IN_T1_WIN_CD
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01
Path_Defs;

symbols = {'ko' 'k*' 'go' 'mo' 'b*' 'r*' 'g*' 'c*'};
lines = {'k-' 'k--' 'g-' 'm-' 'b--' 'r--' 'g--' 'c--'};


%get the column of values of sizes in gratings_params matrix
siz = data.gratings_params(GRAT_WIDTH,:,PATCH1);

%now, get the firing rates for all the trials
spike_rates = data.spike_rates(SpikeChan, :);

%get indices of any NULL conditions (for measuring spontaneous activity)
null_trials = logical( (siz == data.one_time_params(NULL_VALUE)) );

%now, remove trials that do not fall between BegTrial and EndTrial
trials = 1:length(siz);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

%%%% AM added 8/29/16
if exist('skip_trials','var')
    select_trials(skip_trials) = false;
end

%For the size tuning protocol, NULL condition responses are equivalent to
%responses to zero size.  Here we change the NULL_VALUE's to zeros in the
%size vector, so that the zero size point will be included in fit and plot.
%CMA 05/12/05
siz(null_trials) = 0;

stim_type = data.gratings_params(GRAT_TYPE,:,PATCH1);
unique_stim_type = munique(stim_type(~null_trials)');

figure(2);
%set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [450 50 500 573], 'Name', 'Grating Size Tuning Curve');
set(gcf,'PaperPosition', [.51 .51 20 27.5], 'Position', [450 50 500 573], 'Name', 'Grating Size Tuning Curve');
subplot(2, 1, 2);

for j = 1:length(unique_stim_type);

    grating_select = logical((stim_type == unique_stim_type(j)));
    
    plot_x = siz(select_trials & grating_select); %Note: no ~null_trials here because we want the zero size included
    plot_y = spike_rates(select_trials & grating_select);

    %NOTE: Inputs to PlotTuningCurve.m must be column vectors, not row vectors,
    %because of use of munique().
    [px, py, perr] = PlotTuningCurve(plot_x',plot_y',symbols{j},lines{j},1,0);
    
    %Plot the data
    hold on;
    plot(plot_x, plot_y, 'k.')
    errorbar(px, py, perr, perr, symbols{j});
    hold on;
    
    unique_sizes = px; %store category groups for ANOVAs

    %Create inputs to fitting functions
    means = [px py];
    raw = [plot_x' plot_y'];

    %Fit the size tuning curve with an ERROR function (assumes no surround
    %inhibition)
    Erf_pars = erf_fit(means, raw);
    Erf_error = erf_err(Erf_pars);

    %Plot the fitted ERF function
    hold on;
    x_interp = (px(1):0.01:(px(length(px))));
    y_interp = erf_func(x_interp, Erf_pars);
    plot(x_interp, y_interp, 'k-.');
    hold off

    [Erf_chi2(j), Erf_chiP(j)] = Chi2_Test(plot_x, plot_y, 'erf_func', Erf_pars, length(Erf_pars));

    %Fit the size tuning curve with a difference of ERROR functions
    Diff_Erf_pars = diff_erf_fit(means,raw);
    Diff_Erf_error = diff_erf_err(Diff_Erf_pars);

    %Do a sequenctial F test with the ERF and DIFF_ERF fits
    F_SI = ((Erf_error - Diff_Erf_error)/(length(Diff_Erf_pars) - length(Erf_pars)))/(Diff_Erf_error/(length(plot_x)-length(Diff_Erf_pars)));
    P_SI = 1 - fcdf(F_SI, (length(Diff_Erf_pars)-length(Erf_pars)), (length(plot_x)-length(Diff_Erf_pars)));

    %Plot the fitted DIFF_ERF function
    hold on;
    x_interp = (px(1):0.01:(px(length(px))));
    y_interp = diff_erf_func(x_interp,Diff_Erf_pars);
    plot(x_interp, y_interp, 'k-');
    hold off

    %compute discrimination index
    [SDI(j), var_term(j)] = Compute_DDI(plot_x, plot_y);
    
    [Diff_Erf_chi2(j), Diff_Erf_chiP(j)] = Chi2_Test(plot_x, plot_y, 'diff_erf_func', Diff_Erf_pars, length(Diff_Erf_pars));

    [peak_v peak_i] = max(y_interp);

    [val ind] = max(y_interp);
    [min_val min_ind] = min(y_interp(ind:length(y_interp)));

    OptSize(j) = x_interp(peak_i);
    SI(j) = 1-(min_val/val);

    %If P_SI>0.05 (there is no significant surround inhibition), we will set
    %OptSize based on the width parameter of the single erf fit.  SEE GCD,
    %1/30/02
    if (P_SI > 0.05)
        OptSize(j) = 1.163*Erf_pars(2); %this will give us the point that is 90% of the way up the curve
    end

    %now, get the firing rate for NULL condition trials and add spontaneous
    %rate to plot
    null_x = [min(px) max(px)];
    null_rate = mean(data.spike_rates(SpikeChan, null_trials & select_trials));
    null_y = [null_rate null_rate];
    hold on;
    plot(null_x, null_y, 'k--');
    hold off

    y1 = ylim; % 3/14/16 AM corrected from 'YLim'
    ylim([0 y1(2)]); %set the lower limit of the Y axis to zero % % 3/14/16 AM corrected from 'YLim'
    xlabel('Size of Grating (deg)');
    ylabel('Respons (spikes/sec)');

    %Compute R^2 values to characterize goodness of fit
    %Do this using both mean responses and raw (single trial) responses
    R2_pars = Diff_Erf_pars;
    y_fit_mean = diff_erf_func(plot_x', R2_pars);
    y_fit_mean(y_fit_mean < 0) = 0;
    %add a column of ones to yfit to make regress happy
    y_fit_mean = [ones(length(y_fit_mean),1) y_fit_mean];
    [b_mean, bint_mean, r_mean, rint_mean, stats_mean(j,:)] = regress(plot_y', y_fit_mean);

    y_fit_raw = diff_erf_func(plot_x', R2_pars);
    y_fit_raw = [ones(length(y_fit_raw), 1) y_fit_raw];
    [b_raw, bint_raw, r_raw, rint_raw, stats_raw(j,:)] = regress(plot_y', y_fit_raw);
    p_value(j) = spk_anova(plot_y, plot_x, unique_sizes);
    Ke(j) = R2_pars(1);
    a(j) = R2_pars(2);
    Ki(j) = R2_pars(3);
    b_plus_a(j) = R2_pars(4)+R2_pars(2);
    DC(j) = R2_pars(5);
    
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
PrintSizeData(p_value, Ke, a, Ki, b_plus_a, DC, stats_mean, stats_raw, OptSize, SI, SDI, Erf_chi2, Erf_chiP, Diff_Erf_chi2, Diff_Erf_chiP);


output = 1;
if (output == 1)

    %------------------------------------------------------------------------
    %write out all relevant parameters to a cumulative text file, GCD 8/08/01
    %write out one line for each stimu_type for each neuron.
    
    %%%% 5/19/16 AM added conditional below to write into the correct
    %%%% directory if on a specific computer; if on any other lab computer
    %%%% such as 'BURK-LV', the outfile will be written to the usual directory
    if strcmp(getenv('COMPUTERNAME'),'ANDREW') %%% Andrew's MSI laptop name
        outfile = 'C:\Users\AM\Documents\Lab\recordings\GratingSizeSummary.dat';
    else
        outfile = [BASE_PATH 'ProtocolSpecific\GratingSize\GratingSizeSummary.dat'];
    end
    
    printflag = 0;
    if (exist(outfile, 'file') == 0)    %file does not yet exist
        printflag = 1;
    end
    fid = fopen(outfile, 'a');
    if (printflag)
        %%%% 4/15/16 AM edited following line - removed all repeated
        %%%% instances of '\t', which added extra unnecessary columns to
        %%%% output file
        %%% original line: fprintf(fid, 'FILE\t\t\t PrfOR\t PrfSF\t PrfTF\t RFXctr\t RFYctr\t RFDiam\t GrType\t Spont\t SI\t\t AnovaP\t SDI\t\t VarTerm\t MaxX \t MinX \t RsqMn\t PMean\t RsqRaw\t Praw\t\t ErfChi2\t ErfChiP\t DiffErfChi2\t DiffErfChiP\t Ke\t\t a\t\t Ki\t\t\t b+a\t\t DC_off\t Halfmax_slope\t F_SI\t P_SI\t');
        fprintf(fid, 'FILE\t PrfOR\t PrfSF\t PrfTF\t RFXctr\t RFYctr\t RFDiam\t GrType\t Spont\t SI\t AnovaP\t SDI\t VarTerm\t MaxX \t MinX \t RsqMn\t PMean\t RsqRaw\t Praw\t ErfChi2\t ErfChiP\t DiffErfChi2\t DiffErfChiP\t Ke\t a\t Ki\t b+a\t DC_off\t Halfmax_slope\t F_SI\t P_SI\t'); %%%% 4/15/16 AM edited 
        fprintf(fid, '\r\n');
        printflag = 0;
    end
    for j = 1:length(unique_stim_type)
        buff = sprintf('%s\t %5.1f\t %5.2f\t %5.3f\t %5.2f\t %5.2f\t %5.2f\t %5.1f\t %5.4f\t %5.2f\t %5.4f\t %5.2f\t %5.2f\t %5.3f\t %10.8f\t %5.3f\t %6.3f\t %5.2f\t %5.2f\t %7.4f\t %6.4f\t %10.8f\t %6.4f\t %10.8f\t %6.4f\t %10.8f\t %5.2f\t %5.2f\t %6.4f\t %6.4f\t %6.4f\t %6.4f\t %10.8f\t %10.8f\t', ...
            FILE, data.neuron_params(PREFERRED_ORIENTATION, 1), data.neuron_params(PREFERRED_SPATIAL_FREQ, 1), data.neuron_params(PREFERRED_TEMPORAL_FREQ, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
            unique_stim_type(j), null_rate, SI(j), p_value(j), SDI(j), var_term(j), max(unique_sizes), min(unique_sizes), stats_mean(j,1), stats_mean(j,3), stats_raw(j,1), stats_raw(j,3), Erf_chi2(j), Erf_chiP(j), Diff_Erf_chi2(j), Diff_Erf_chiP(j), Ke(j), a(j), Ki(j), b_plus_a(j), DC(j), halfmax_slope(j), F_SI, P_SI );
        fprintf(fid, '%s', buff);
        fprintf(fid, '\r\n');
    end
    fclose(fid);
    %------------------------------------------------------------------------

end  %if (output == 1)


return;