%-----------------------------------------------------------------------------------------------------------------------
%-- TFTuningCurve.m -- Plots a TF tuning curve,for sine
%and square waves
%--	GCD, 1/12/06, CMA 4/06
%-----------------------------------------------------------------------------------------------------------------------

%%%% 4/27/16 AM added conditional when writing outfile to write into the correct
%%%% directory if on a specific computer; functionality should be
%%%% unchanged if using on any lab computer
%%%% 4/27/16 AM edited formatting of output file column
%%% 3/29/16 AM added escape from endless while loop..... endless loops
%%%%    seem to occur when the cell is not well-tuned for the parameter
%%%%    (function does not fit the data well)
% 3/14/16 AM fixed typo
%%%% 8/23/16 AM added optional 'skip_trials' argin - these trials will not be analyzed;
%%%%    trials not included within BegTrial-EndTrial will still be skipped

function TFTuningCurve(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, skip_trials)

countlimiter = 1; %%% AM added 3/29/16 - toggle on to limit while loops to countlimit cycles
countlimit = 1e6; %%% AM added 3/29/16 to break endless loop

TEMPO_Defs;		%needed for defines like IN_T1_WIN_CD
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01
Path_Defs;

symbols = {'ko' 'k*' 'go' 'mo' 'b*' 'r*' 'g*' 'c*'};
lines = {'k-' 'k--' 'g-' 'm-' 'b--' 'r--' 'g--' 'c--'};


%get the column of values of SFs in gratings_params matrix
tf = data.gratings_params(GRAT_TEMPORAL_FREQ,:,PATCH1);

%get indices of any NULL conditions (for measuring spontaneous activity)
null_trials = logical( (tf == data.one_time_params(NULL_VALUE)) );

%get the column of stimulus type values
stim_type = data.gratings_params(GRAT_TYPE,:,PATCH1);
unique_stim_type = munique(stim_type(~null_trials)');

%now, get the firing rates for all the trials
spike_rates = data.spike_rates(SpikeChan, :);

unique_tf = munique(tf(~null_trials)');

%now, remove trials that do not fall between BegTrial and EndTrial
trials = 1:length(tf);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

%%%% AM added 8/29/16
if exist('skip_trials','var')
    select_trials(skip_trials) = false;
end

figure(2);
set(gcf,'PaperPosition', [.51 .51 20 27.5], 'Position', [450 50 500 573], 'Name', 'Grating Temporal Frequency Tuning Curve');
subplot(2, 1, 2);

%Code below this point edited by CMA 04/17/06

for j = 1:length(unique_stim_type)

    grating_select = logical( (stim_type == unique_stim_type(j)) );
    
    plot_x=tf(~null_trials & select_trials & grating_select);
    plot_y=spike_rates(~null_trials & select_trials & grating_select);

    %NOTE:  Inputs to PlotTuningCurve.m must be column vectors, not row vectors
    %because of munique().
    [px, py, perr] = PlotTuningCurve(plot_x', plot_y', symbols{j}, lines{j}, 1, 0);
    %NOTE: last argument=0 means just get output, no plot

    unique_tfreqs=px; %store category groups for ANOVAs

    %Plot the data
    hold on;
    plot(plot_x, plot_y, 'k.');
    errorbar(px, py, perr, perr, symbols{j});
    hold on;
    if (j==1)
        sine_px=px; sine_py=py; sine_perr=perr;
    end
    if (j==2)
        sq_px=px; sq_py=py; sq_perr=perr;
    end

    %compute discrimination index
    [TFDI(j), var_term(j)] = Compute_DDI(plot_x, plot_y);

    %Now, fit the data with a log gaussian curve and plot this as well
    means = [px py];
    raw = [plot_x' plot_y'];
    [pars] = TF_loggaussfit(means, raw) %last arg=0 => only allow positive going fit
    x_interp = (px(1)):0.01:(px(length(px)));
    y_interp = TF_loggaussfunc(x_interp, pars);
    plot(x_interp, y_interp, 'k-');

    if (j==1)
        sine_x_interp=x_interp; sine_y_interp=y_interp;
    end
    if (j==2)
        sq_x_interp=x_interp; sq_y_interp=y_interp;
    end

    %now, get the firing rate for NULL condition trials and add spontaneous
    %rate to plot
    null_x = [min(px) max(px)];
    null_resp = data.spike_rates(SpikeChan, null_trials & select_trials);
    null_rate = mean(null_resp);
    null_y = [null_rate null_rate];
    hold on;
    plot(null_x, null_y, 'k--');
    hold off;

    if (j==1)
        sine_null_x=null_x; sine_null_y=null_y;
    end
    if (j==2)
        sq_null_x=null_x; sq_null_y=null_y;
    end
    
    yl=ylim;
    ylim([0 yl(2)]); %set lower limit of the Y axis to zero
    xlabel('Grating Temporal Frequency (cyc/sec)');  % 3/14/16 AM changed from '(cyc/deg)'
    ylabel('Response (spikes/sec)');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Added on 04/21/06, suppressed until TF_loggaussfit is working
    %Compute R^2 for both means and raw values
    y_fit = TF_loggaussfunc(px, pars);
    y_fit(y_fit<0) = 0;
    %add a column of ones to yfit to make regress happy
    y_fit = [ones(length(y_fit), 1) y_fit];
    [b, bint, r, rint, stats1(j,:)] = regress(py, y_fit);

    y_fit_raw = TF_loggaussfunc(plot_x', pars);
    y_fit_raw(y_fit_raw<0) = 0;
    y_fit_raw = [ones(length(y_fit_raw), 1) y_fit_raw];
    [b, bint, r, rint, stats2(j,:)] = regress(plot_y', y_fit_raw);

    %Calculate some metrics and stats to print in the plot

    p_value(j) = spk_anova(plot_y, plot_x, unique_tfreqs);
    base_rate(j) = pars(1);
    amplitude(j) = pars(2);
    peak(j) = pars(3);
    st_dev(j) = pars(4);
    log_offset(j) = pars(5);
    max_rate(j) = base_rate(j) + amplitude(j);

    TFMI(j) = 1 - (base_rate - null_rate)/(max_rate - null_rate);

    %Do chi-square goodness of fit test
    [chi2(j), chiP(j)] = Chi2_Test(plot_x, plot_y, 'TF_loggaussfunc', pars, length(pars))

    %Rsq_means = stats1(1);
    %P_means = stats1(3);
    %Rsq_raw = stats2(1);
    %P_raw = stats2(3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%% AM note: get right sided HWHM (search from response peak to
        %%% right side of interpolated values in tuning curve
    [pmax{j}.y max_i] = max(y_interp);
    pmax{j}.x = x_interp(max_i);
    [pmin{j}.y min_i] = min(y_interp);
    pmin{j}.x = x_interp(min_i);

    index = max_i;
    ysearch = pmax{j}.y;
    count = 0; %%% AM added 3/29/16 to break endless loop
    while (ysearch > (pmin{j}.y + (pmax{j}.y-pmin{j}.y)/2))
        count = count+1; %%% AM added 3/29/16 to break endless loop 
        if countlimiter && count > countlimit; %%% AM added 3/29/16 to break endless loop
            break  %%% AM added 3/29/16 to break endless loop
        end   %%% AM added 3/29/16 to break endless loop
        ysearch = y_interp(index); 
        if (index == length(y_interp))
            break;
        else
            index = index + 1;
        end
    end
    tf_width_halfmax(j) = x_interp(index) - pmax{j}.x;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find widths at various percentages of the difference between the max
    % and min.  Send these to an output file so that they can be compared.
    index = max_i;
    ysearch = pmax{j}.y;
    count = 0; %%% AM added 3/29/16 to break endless loop
    while (ysearch > (pmin{j}.y + (pmax{j}.y-pmin{j}.y)*0.25))
        count = count+1; %%% AM added 3/29/16 to break endless loop 
        if countlimiter && count > countlimit; %%% AM added 3/29/16 to break endless loop
            break  %%% AM added 3/29/16 to break endless loop
        end   %%% AM added 3/29/16 to break endless loop
        ysearch = y_interp(index);
        if (index == length(y_interp))
            tf_width_25_percent(j) = 9999;
        else
            index = index + 1;
        end
    end
    tf_width_25_percent(j) = x_interp(index) - pmax{j}.x;
    
    index = max_i;
    ysearch = pmax{j}.y;
    count = 0; %%% AM added 3/29/16 to break endless loop
    while (ysearch > (pmin{j}.y + (pmax{j}.y-pmin{j}.y)*0.20))
        count = count+1; %%% AM added 3/29/16 to break endless loop 
        if countlimiter && count > countlimit; %%% AM added 3/29/16 to break endless loop
            break  %%% AM added 3/29/16 to break endless loop
        end   %%% AM added 3/29/16 to break endless loop
        ysearch = y_interp(index);
        if (index == length(y_interp))
            tf_width_20_percent(j) = 9999;
        else
            index = index + 1;
        end
    end
    tf_width_20_percent(j) = x_interp(index) - pmax{j}.x;
    
    index = max_i;
    ysearch = pmax{j}.y;
    count = 0; %%% AM added 3/29/16 to break endless loop 
    while (ysearch > (pmin{j}.y + (pmax{j}.y-pmin{j}.y)*0.15))
        count = count+1; %%% AM added 3/29/16 to break endless loop 
        if countlimiter && count > countlimit; %%% AM added 3/29/16 to break endless loop
            break  %%% AM added 3/29/16 to break endless loop
        end   %%% AM added 3/29/16 to break endless loop
        ysearch = y_interp(index);
        if (index == length(y_interp))
            tf_width_15_percent(j) = 9999;
        else
            index = index + 1;
        end
    end
    tf_width_15_percent(j) = x_interp(index) - pmax{j}.x;
    
    index = max_i;
    ysearch = pmax{j}.y;
    count = 0; 
    while (ysearch > (pmin{j}.y + (pmax{j}.y-pmin{j}.y)*0.10))
        count = count+1; %%% AM added 3/29/16 to break endless loop 
        if countlimiter && count > countlimit; %%% AM added 3/29/16 to break endless loop
            break  %%% AM added 3/29/16 to break endless loop
        end   %%% AM added 3/29/16 to break endless loop
        ysearch = y_interp(index);
        if (index == length(y_interp))
            tf_width_10_percent(j) = 9999;
        else
            index = index + 1;
        end
    end
    tf_width_10_percent(j) = x_interp(index) - pmax{j}.x;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

%now, print out some useful information in the upper subplot
subplot(2, 1, 1);
PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

%now print out useful values
% pmax, pmin, py,
PrintTemporalFreqData(p_value, base_rate, amplitude, peak, st_dev, log_offset, max_rate, stats1, stats2, chi2, chiP, TFDI, tf_width_halfmax);

output = 1;
if (output == 1)
% 
%     i = size(PATH,2) - 1;
%     while PATH(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
%         i = i - 1;
%     end   
%     PATHOUT = [PATH(1:i) 'Analysis\Tuning\'];
%     i = size(FILE,2) - 1;
%     while FILE(i) ~='.'
%         i = i - 1;
%     end
% %also write out data in form suitable for plotting tuning curve with Origin.
%     FILEOUT2 = [FILE(1:i) 'grat_tf_curv'];
%     fileid = [PATHOUT FILEOUT2];
%     proffid = fopen(fileid, 'w');
%     fprintf(proffid,'TFIn\tFit\tTF\tAvgResp\tStdErr\tTF2\tSpon\n');
%     for i=1:length(x_interp)
%         fprintf(proffid,'%6.4f\t%6.4f\t', sine_x_interp(i), sine_y_interp(i));
%         if (i <= length(px))
%             fprintf(proffid,'%6.4f\t%6.2f\t%6.3f\t', sine_px(i), sine_py(i), sine_perr(i));
%         else
%             fprintf(proffid,'\t\t\t');
%         end
%         if (i <= 2)
%             fprintf(proffid,'%6.4f\t%6.2f\n',sine_null_x(i),sine_null_y(i));
%         else
%             fprintf(proffid,'\t\n');
%         end
%     end
%     fclose(proffid);
    
    %------------------------------------------------------------------------
    %write out all relevant parameters to a cumulative text file, GCD 8/08/01
    %write out one line for each stimu_type for each neuron.
    
    %%%% 4/27/16 AM added conditional below to write into the correct
    %%%% directory if on a specific computer; if on any other lab computer
    %%%% such as 'BURK-LV', the outfile will be written to the usual directory
    if strcmp(getenv('COMPUTERNAME'),'ANDREW') %%% Andrew's MSI laptop name
        outfile = 'C:\Users\AM\Documents\Lab\recordings\GratingTF_Summary.dat';
    else
        outfile = [BASE_PATH 'ProtocolSpecific\GratingTemporalFreq\GratingTF_Summary.dat'];
    end
    
    printflag = 0;
    if (exist(outfile, 'file') == 0)    %file does not yet exist
        printflag = 1;
    end
    fid = fopen(outfile, 'a');
    if (printflag)
        %%%% 4/27/16 AM edited following line - removed all repeated
        %%%% instances of '\t', which added extra unnecessary columns to
        %%%% output file
        %%% original line: fprintf(fid, 'FILE\t\t\t PrfOR\t PrfSF\t PrfTF\t RFXctr\t RFYctr\t RFDiam\t GrType\t MaxLoc\t MaxRsp\t MinLoc\t MinRsp\t Spont\t TFMI\t\t AnovaP\t TFDI  \t VarTerm\t MaxX \t MinX \t TFHWHM\t TFHW25\t TFHW20\t TFHW15\t TFHW10\t RsqMn\t PMean\t RsqRaw\t Praw\t\t Chi2\t\t ChiP\t\t Brat\t\t Ampl\t\t Peak\t\t SD\t\t Offst\t\t');
        fprintf(fid, 'FILE\t PrfOR\t PrfSF\t PrfTF\t RFXctr\t RFYctr\t RFDiam\t GrType\t MaxLoc\t MaxRsp\t MinLoc\t MinRsp\t Spont\t TFMI\t AnovaP\t TFDI  \t VarTerm\t MaxX \t MinX \t TFHWHM\t TFHW25\t TFHW20\t TFHW15\t TFHW10\t RsqMn\t PMean\t RsqRaw\t Praw\t Chi2\t ChiP\t Brat\t Ampl\t Peak\t SD\t Offst\t\t');
        fprintf(fid, '\r\n');
        printflag = 0;
    end
    for j = 1:length(unique_stim_type)
        buff = sprintf('%s\t %5.1f\t %5.2f\t %5.3f\t %5.2f\t %5.2f\t %5.2f\t %5.1f\t %5.4f\t %5.2f\t %5.4f\t %5.2f\t %5.2f\t %5.3f\t %10.8f\t %5.3f\t %6.3f\t %5.2f\t %5.2f\t %7.4f\t %7.4f\t %7.4f\t %7.4f\t %7.4f\t %6.4f\t %10.8f\t %6.4f\t %10.8f\t %6.4f\t %10.8f\t %5.2f\t %5.2f\t %6.4f\t %6.4f\t %6.4f\t', ...
            FILE, data.neuron_params(PREFERRED_ORIENTATION, 1), data.neuron_params(PREFERRED_SPATIAL_FREQ, 1), data.neuron_params(PREFERRED_TEMPORAL_FREQ, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
            unique_stim_type(j), pmax{j}.x, pmax{j}.y, pmin{j}.x, pmin{j}.y, null_rate, TFMI(j), p_value(j), TFDI(j), var_term(j), max(unique_tf), min(unique_tf), tf_width_halfmax(j), tf_width_25_percent(j), tf_width_20_percent(j), tf_width_15_percent(j), tf_width_10_percent(j), stats1(j,1), stats1(j,3), stats2(j,1), stats2(j,3), chi2(j), chiP(j), base_rate(j), amplitude(j), peak(j), st_dev(j), log_offset(j) );
        fprintf(fid, '%s', buff);
        fprintf(fid, '\r\n');
    end
    fclose(fid);
    %------------------------------------------------------------------------

end  %if (output == 1)

return;