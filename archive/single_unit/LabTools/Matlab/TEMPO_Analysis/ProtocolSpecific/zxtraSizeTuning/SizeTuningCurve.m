%-----------------------------------------------------------------------------------------------------------------------
%-- SizeTuningCurve.m -- Plots a size tuning curve and computes/plots a difference of gaussian fit
%--	JDN, 2/18/00; modified extensively by GCD, 11/2001
%-----------------------------------------------------------------------------------------------------------------------
function SizeTuningCurve(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;
Path_Defs;

%get the column of values of directions in the dots_params matrix
size_vals = data.dots_params(DOTS_AP_XSIZ,:,PATCH1);

%now, get the firing rates for all the trials 
spike_rates = data.spike_rates(SpikeChan, :);

%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (size_vals == data.one_time_params(NULL_VALUE)) );

%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(size_vals);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

stim_type = data.dots_params(GRAT_TYPE,:,PATCH1);
unique_stim_type = munique(stim_type(~null_trials)');

%for the size tuning protocol, NULL condition responses are equivalent to responses to zero size.
%here we change the NULL_VALUE's to zeros in the size vector, so that the zero size point will be 
%included in fit and plot.  GCD, 10/16/01
size_vals(null_trials) = 0;

unique_size_vals = munique(size_vals');

figure;
set(gcf,'PaperPosition', [.51 .51 20 27.5], 'Position', [250 250 500 573], 'Name', 'Size Tuning Curve');
subplot(2, 1, 2);

for j = 1:length(unique_stim_type)

    grating_select = logical((stim_type == unique_stim_type(j)));
    
    plot_x = size_vals(select_trials & grating_select);
    plot_y = spike_rates(select_trials & grating_select);

    %NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
    [px, py, perr] = PlotTuningCurve(plot_x', plot_y', symbols{j}, lines{j}, 1, 0);
    %NOTE: last arg=0 means just get output, no plot
    
    unique_sizes = px; %store category group for ANOVAs
    
    hold on;
    plot(plot_x, plot_y, 'k.');
    errorbar(px, py, perr, perr, symbols{j});
    hold on;

    %compute discrimination index
    [SDI(j), var_term(j)] = Compute_DDI(plot_x, plot_y);

    %create inputs to fitting functions
    means = [px py];
    raw = [plot_x' plot_y'];

    %Fit the size tuning curve with an ERROR function (assumes no surround inhibition)
    Erf_pars = erf_fit(means,raw);
    Erf_error = erf_err(Erf_pars);

    %Plot the fitted ERF function
    hold on;
    x_interp = (px(1):.01:(px(length(px))));
    y_interp = erf_func(x_interp,Erf_pars);
    plot(x_interp, y_interp, 'k.');
    hold off

    [Erf_chi2(j), Erf_chiP(j)] = Chi2_Test(plot_x, plot_y, 'erf_func', Erf_pars, length(Erf_pars))

    %Fit the size tuning curve with a difference of ERROR functions.
    Diff_Erf_pars = diff_erf_fit(means,raw);
    Diff_Erf_error = diff_erf_err(Diff_Erf_pars);

    %Do a sequential F test with the ERF and DIFF_ERF fits
    F_SI = ( (Erf_error - Diff_Erf_error)/(length(Diff_Erf_pars)-length(Erf_pars)) ) / ( Diff_Erf_error/(length(plot_x)-length(Diff_Erf_pars)) );
    P_SI = 1 - fcdf(F_SI, (length(Diff_Erf_pars)-length(Erf_pars)), (length(plot_x)-length(Diff_Erf_pars)) );

    %Plot the fitted DIFF_ERF function
    hold on;
    x_interp = (px(1):.01:(px(length(px))));
    y_interp = diff_erf_func(x_interp,Diff_Erf_pars);
    plot(x_interp, y_interp, 'k-');
    hold off

    [Diff_Erf_chi2(j), Diff_Erf_chiP(j)] = Chi2_Test(plot_x, plot_y, 'diff_erf_func', Diff_Erf_pars, length(Diff_Erf_pars));

    [peak_v peak_i] = max(y_interp);

    [val ind] = max(y_interp);
    [min_val min_ind] = min(y_interp(ind:length(y_interp)));

    OptSize(j) = x_interp(peak_i);
    SI(j) = 1-(min_val/val);

    %If P_SI>0.05 (there is no significant surround inhibition), we will set OptSize based on the
    %width parameter of the single erf fit.  GCD, 1/30/02
    if (P_SI > 0.05)
        OptSize(j) = 1.163*Erf_pars(2);  %this will give us the point that is 90% of the way up the curve
    end

    %now, get the firing rate for NULL condition trials and add spontaneous rate to plot
    null_x = [min(px) max(px)];
    null_rate = mean(data.spike_rates(SpikeChan, null_trials & select_trials));
    null_y = [null_rate null_rate];
    hold on;
    plot(null_x, null_y, 'k--');
    hold off;

    yl = YLim;
    YLim([0 yl(2)]);	% set the lower limit of the Y axis to zero
    XLabel('Size of Aperture (deg)');
    YLabel('Response (spikes/sec)');

    %Compute R^2 values to characterize goodness of fit
    %Do this using both the mean responses and raw (single trial) responses
    R2_pars = Diff_Erf_pars;
    y_fit_mean = diff_erf_func(px, R2_pars);
    y_fit_mean(y_fit_mean < 0) = 0;
    %add a column of ones to yfit to make regress happy
    y_fit_mean = [ones(length(y_fit_mean),1) y_fit_mean];
    [b_mean(j), bint_mean(j), r_mean(j), rint_mean(j), stats_mean(j)] = regress(py, y_fit_mean);

    y_fit_raw = diff_erf_func(plot_x', R2_pars);
    y_fit_raw(y_fit_raw < 0) = 0;
    y_fit_raw = [ones(length(y_fit_raw),1) y_fit_raw];
    [b_raw(j), bint_raw(j), r_raw(j), rint_raw(j), stats_raw(j)] = regress(plot_y', y_fit_raw);
    
    p_value(j) = spk_anova(plot_y, plot_x, unique_sizes);

end


% %now, print out some useful information in the upper subplot
% subplot(2, 1, 1);
% PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
% 
% %now print out useful values
% % pmax, pmin, py,
% PrintSizeData(p_value, base_rate, amplitude, peak, st_dev, log_offset, max_rate, stats1, stats2, chi2, chiP, SFDI, sf_width_halfmax);
% 
% output = 1;
% if (output == 1)
% 
%     %------------------------------------------------------------------------
%     %write out all relevant parameters to a cumulative text file, GCD 8/08/01
%     %write out one line for each stimu_type for each neuron.
%     outfile = [BASE_PATH 'ProtocolSpecific\GratingSpatialFreq\GratingSF_Summary.dat'];
%     printflag = 0;
%     if (exist(outfile, 'file') == 0)    %file does not yet exist
%         printflag = 1;
%     end
%     fid = fopen(outfile, 'a');
%     if (printflag)
%         fprintf(fid, 'FILE\t\t\t PrfOR\t PrfSF\t PrfTF\t RFXctr\t RFYctr\t RFDiam\t GrType\t MaxLoc\t MaxRsp\t MinLoc\t MinRsp\t Spont\t SFMI\t\t AnovaP\t SFDI  \t VarTerm\t MaxX \t MinX \t SFHWHM\t RsqMn\t PMean\t RsqRaw\t Praw\t\t Chi2\t\t ChiP\t\t Brat\t\t Ampl\t\t Peak\t\t SD\t\t Offst\t\t');
%         fprintf(fid, '\r\n');
%         printflag = 0;
%     end
%     for j = 1:length(unique_stim_type)
%         buff = sprintf('%s\t %5.1f\t %5.2f\t %5.3f\t %5.2f\t %5.2f\t %5.2f\t %5.1f\t %5.4f\t %5.2f\t %5.4f\t %5.2f\t %5.2f\t %5.3f\t %10.8f\t %5.3f\t %6.3f\t %5.2f\t %5.2f\t %7.4f\t %6.4f\t %10.8f\t %6.4f\t %10.8f\t %6.4f\t %10.8f\t %5.2f\t %5.2f\t %6.4f\t %6.4f\t %6.4f\t', ...
%             FILE, data.neuron_params(PREFERRED_ORIENTATION, 1), data.neuron_params(PREFERRED_SPATIAL_FREQ, 1), data.neuron_params(PREFERRED_TEMPORAL_FREQ, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
%             unique_stim_type(j), pmax{j}.x, pmax{j}.y, pmin{j}.x, pmin{j}.y, null_rate, SFMI(j), p_value(j), SFDI(j), var_term(j), max(unique_sf), min(unique_sf), sf_width_halfmax(j), stats1(j,1), stats1(j,3), stats2(j,1), stats2(j,3), chi2(j), chiP(j), base_rate(j), amplitude(j), peak(j), st_dev(j), log_offset(j) );
%         fprintf(fid, '%s', buff);
%         fprintf(fid, '\r\n');
%     end
%     fclose(fid);
%     %------------------------------------------------------------------------
% 
% end  %if (output == 1)


return;


% %now, print out some useful information in the upper subplot
% subplot(2, 1, 1);
% PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
% subplot(2, 1, 1);
% diff_erf_parsvalues = sprintf('Ke = %4.4f, a = %4.4f\nKi = %4.4f, b = %4.4f\nR0 = %5.3f\nOptSize = %4.1f, SI = %4.3f pct\nFseq = %6.4f, Pseq = %8.6f\nRsqraw = %6.4f  Praw = %10.8f\nRsqmean = %6.4f Pmean = %10.8f', Diff_Erf_pars, OptSize, SI*100, F_SI, P_SI, stats_raw(1), stats_raw(3), stats_mean(1), stats_mean(3));
% axval = axis;
% text(-10, axval(3) +2, diff_erf_parsvalues, 'FontSize', 8);
% erf_parsvalues = sprintf('K = %4.4f\na = %4.4f\nR0 = %5.3f\n', Erf_pars);
% text(70, axval(3) +4, erf_parsvalues, 'FontSize', 8);
% 
% 
% 
% 
% %------------------------------------------------------------------------
% %write out all relevant parameters to a cumulative text file, GCD 10/18/01
% outfile = [BASE_PATH 'ProtocolSpecific\SizeTuning\SizeTuningSummary.dat'];
% %PATHOUT = 'Z:\Users\Jerry\GradAnalysis\';
% %outfile = [PATHOUT 'GradSI_metrics_013003.dat'];
% printflag = 0;
% if (exist(outfile, 'file') == 0)    %file does not yet exist
%     printflag = 1;
% end
% fid = fopen(outfile, 'a');
% if (printflag)
%     fprintf(fid, 'FILE\t\t PrDir\t PrSpd\t PrHDsp\t RFX\t RFY\t RFDiam\t K\t a\t R0\t Ke\t ae\t Ki\t bi\t R0\t OptSiz\t PctSI\t Fseq\t Pseq\t\t R2raw\t Praw\t\t R2mean\t Pmean\t\t Chi2E\t\t ChiPE\t\t Chi2DE\t\t ChiPDE\t\t');
%     fprintf(fid, '\r\n');
%     printflag = 0;
% end
% buff = sprintf('%s\t %6.1f\t %6.2f\t %6.3f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %5.3f\t %5.3f\t %6.2f\t %5.3f\t %6.2f\t%5.3f\t %5.3f\t %6.2f\t %6.2f\t %6.2f\t %10.8f\t %6.4f\t %10.8f\t %6.4f\t %10.8f\t %8.5f\t %10.8f\t %8.5f\t %10.8f\t', ...
%     FILE, data.neuron_params(PREFERRED_DIRECTION, 1), data.neuron_params(PREFERRED_SPEED, 1), data.neuron_params(PREFERRED_HDISP, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
%     Erf_pars, Diff_Erf_pars, OptSize, SI*100, F_SI, P_SI, stats_raw(1), stats_raw(3), stats_mean(1), stats_mean(3), Erf_chi2, Erf_chiP, Diff_Erf_chi2, Diff_Erf_chiP);
% fprintf(fid, '%s', buff);
% fprintf(fid, '\r\n');
% fclose(fid);
% %------------------------------------------------------------------------




% %------------------------------------------------------------------------
% %write out data in form suitable for plotting tuning curve with Origin.
% i = size(PATH,2) - 1;
% while PATH(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
%     i = i - 1;
% end   
% PATHOUT = [PATH(1:i) 'Analysis\Tuning\'];
% i = size(FILE,2) - 1;
% while FILE(i) ~='.'
%     i = i - 1;
% end
% 
% FILEOUT2 = [FILE(1:i) 'size_curv_fit'];
% fileid = [PATHOUT FILEOUT2];
% proffid = fopen(fileid, 'w');
% fprintf(proffid,'SizIn\tFit\tSize\tAvgResp\tStdErr\tSiz2\tSpon\n');
% for i=1:length(x_interp)
%     fprintf(proffid,'%6.2f\t%6.2f\t', x_interp(i), y_interp(i));
%     if (i <= length(px))
%         fprintf(proffid,'%6.2f\t%6.2f\t%6.3f\t', px(i), py(i), perr(i));
%     else
%         fprintf(proffid,'\t\t\t');
%     end
%     if (i <= 2)
%         fprintf(proffid,'%6.2f\t%6.2f\n',null_x(i),null_y(i));
%     else
%         fprintf(proffid,'\t\n');
%     end
% end
% fclose(proffid);


% grad_print = 0;
% if (grad_print == 1)
%     PATHOUT = 'Z:\Users\Jerry\GradAnalysis\';
%     
%     line = sprintf('%s %1.12f', FILE, SI);	
%     outfile = [PATHOUT 'SI_metric.dat'];
%     printflag = 0;
%     if (exist(outfile, 'file') == 0)    %file does not yet exist
%         printflag = 1;
%     end
%     fid = fopen(outfile, 'a');
%     if (printflag)
%         fprintf(fid, 'FILE         SI');
%         fprintf(fid, '\r\n');
%     end
%     fprintf(fid, '%s', [line]);
%     fprintf(fid, '\r\n');
%     fclose(fid);
% end

