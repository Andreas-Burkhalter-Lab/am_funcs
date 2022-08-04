%Edited on 10/12/07 by RLS to plot against log-SF
%-----------------------------------------------------------------------------------------------------------------------
%-- SFTuningCurve.m -- Plots an SF tuning curve,for sine
%and square waves
%--	GCD, 1/12/06
%-----------------------------------------------------------------------------------------------------------------------
function SFTuningCurve(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;		%needed for defines like IN_T1_WIN_CD
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01
Path_Defs;

symbols = {'ko' 'k*' 'go' 'mo' 'b*' 'r*' 'g*' 'c*'};
lines = {'k-' 'k--' 'g-' 'm-' 'b--' 'r--' 'g--' 'c--'};


if length(intersect(unique(data.laser_values),[0 1]))==1
    %% DHK - 2 Laser Conditon -> one experiments
    %get the column of values of SFs in gratings_params matrix
    sf = data.gratings_params(GRAT_SPATIAL_FREQ,:,PATCH1);
    
    %get indices of any NULL conditions (for measuring spontaneous activity)
    null_trials = logical( (sf == data.one_time_params(NULL_VALUE)) );
    
    %get the column of stimulus type values
    stim_type = data.gratings_params(GRAT_TYPE,:,PATCH1);
    unique_stim_type = munique(stim_type(~null_trials)');
    
    %now, get the firing rates for all the trials
    spike_rates = data.spike_rates(SpikeChan, :);
    
    unique_sf = munique(sf(~null_trials)')
    
    %now, remove trials that do not fall between BegTrial and EndTrial
    trials = 1:length(sf);		% a vector of trial indices
    select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );
    
    figure(2);
    set(gcf,'PaperPosition', [.51 .51 20 27.5], 'Position', [450 50 500 573], 'Name', 'Grating Spatial Frequency Tuning Curve');
    subplot(2, 1, 2);
    
    %Code below this point edited by CMA 04/17/06
    
    for j = 1:length(unique_stim_type);
        
        grating_select = logical( (stim_type == unique_stim_type(j)) );
        
        plot_x=sf(~null_trials & select_trials & grating_select);
        plot_y=spike_rates(~null_trials & select_trials & grating_select);
        
        %NOTE:  Inputs to PlotTuningCurve.m must be column vectors, not row vectors
        %because of munique().
        [px, py, perr] = PlotTuningCurve(plot_x', plot_y', symbols{j}, lines{j}, 1, 0);
        %NOTE: last argument=0 means just get output, no plot
        
        unique_sfreqs=px; %store category groups for ANOVAs
        
        %Plot the data
        hold on;
        plot(plot_x, plot_y, 'k.');
        errorbar(px, py, perr, perr, symbols{j});
        hold on;
        
        %compute discrimination index
        [SFDI(j), var_term(j)] = Compute_DDI(plot_x, plot_y);
        
        %Now, fit the data with a Gaussian curve and plot this as well
        means = [px py];
        raw = [plot_x' plot_y'];
        [pars] = SF_loggaussfit(means, raw); %last arg=0 => only allow positive going fit
        x_interp = (px(1)):0.001:(px(length(px)));
        y_interp = SF_loggaussfunc(x_interp, pars);
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
        xlabel('Grating Spatial Frequency (cyc/deg)');
        ylabel('Response (spikes/sec)');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Added on 04/21/06, suppressed until SF_loggaussfit is working
        %Compute R^2 for both means and raw values
        y_fit = SF_loggaussfunc(px, pars);
        y_fit(y_fit<0) = 0;
        %add a column of ones to yfit to make regress happy
        y_fit = [ones(length(y_fit), 1) y_fit];
        [b, bint, r, rint, stats1(j,:)] = regress(py, y_fit);
        
        y_fit_raw = SF_loggaussfunc(plot_x', pars);
        y_fit_raw(y_fit_raw<0) = 0;
        y_fit_raw = [ones(length(y_fit_raw), 1) y_fit_raw];
        [b, bint, r, rint, stats2(j,:)] = regress(plot_y', y_fit_raw);
        
        %Calculate some metrics and stats to print in the plot
        
        p_value(j) = spk_anova(plot_y, plot_x, unique_sfreqs);
        base_rate(j) = pars(1);
        amplitude(j) = pars(2);
        peak(j) = pars(3);
        st_dev(j) = pars(4);
        log_offset(j) = pars(5);
        max_rate(j) = base_rate(j) + amplitude(j);
        
        SFMI(j) = 1 - (base_rate - null_rate)/(max_rate - null_rate);
        
        %Do chi-square goodness of fit test
        [chi2(j), chiP(j)] = Chi2_Test(plot_x, plot_y, 'SF_loggaussfunc', pars, length(pars));
        
        %Rsq_means = stats1(1);
        %P_means = stats1(3);
        %Rsq_raw = stats2(1);
        %P_raw = stats2(3);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [pmax{j}.y max_i] = max(y_interp);
        pmax{j}.x = x_interp(max_i);
        [pmin{j}.y min_i] = min(y_interp);
        pmin{j}.x = x_interp(min_i);
        
        index = max_i;
        ysearch = pmax{j}.y;
        while (ysearch > (pmin{j}.y + (pmax{j}.y-pmin{j}.y)/2))
            ysearch = y_interp(index);
            if (index == length(y_interp))
                break;
            else
                index = index + 1;
            end
        end
        sf_width_halfmax(j) = x_interp(index) - pmax{j}.x;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Find widths at various percentages of the difference between the max
        % and min.  Send these to an output file so that they can be compared.
        index = max_i;
        ysearch = pmax{j}.y;
        while (ysearch > (pmin{j}.y + (pmax{j}.y-pmin{j}.y)*0.25))
            ysearch = y_interp(index);
            if (index == length(y_interp))
                sf_width_25_percent(j) = 9999;
            else
                index = index + 1;
            end
        end
        sf_width_25_percent(j) = x_interp(index) - pmax{j}.x;
        
        index = max_i;
        ysearch = pmax{j}.y;
        while (ysearch > (pmin{j}.y + (pmax{j}.y-pmin{j}.y)*0.20))
            ysearch = y_interp(index);
            if (index == length(y_interp))
                sf_width_20_percent(j) = 9999;
            else
                index = index + 1;
            end
        end
        sf_width_20_percent(j) = x_interp(index) - pmax{j}.x;
        
        index = max_i;
        ysearch = pmax{j}.y;
        while (ysearch > (pmin{j}.y + (pmax{j}.y-pmin{j}.y)*0.15))
            ysearch = y_interp(index);
            if (index == length(y_interp))
                sf_width_15_percent(j) = 9999;
            else
                index = index + 1;
            end
        end
        sf_width_15_percent(j) = x_interp(index) - pmax{j}.x;
        
        index = max_i;
        ysearch = pmax{j}.y;
        while (ysearch > (pmin{j}.y + (pmax{j}.y-pmin{j}.y)*0.10))
            ysearch = y_interp(index);
            if (index == length(y_interp))
                sf_width_10_percent(j) = 9999;
            else
                index = index + 1;
            end
        end
        sf_width_10_percent(j) = x_interp(index) - pmax{j}.x;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
    end
    
    %now, print out some useful information in the upper subplot
    subplot(2, 1, 1);
    PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    
    %now print out useful values
    % pmax, pmin, py,
    PrintSpatialFreqData(p_value, base_rate, amplitude, peak, st_dev, log_offset, max_rate, stats1, stats2, chi2, chiP, SFDI, sf_width_halfmax);
    
    output = 1;
    if (output == 1)
        
        i = size(PATH,2) - 1;
        while PATH(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
            i = i - 1;
        end
        PATHOUT = [PATH(1:i) 'Analysis\Tuning\'];
        i = size(FILE,2) - 1;
        while FILE(i) ~='.'
            i = i - 1;
        end
        %also write out data in form suitable for plotting tuning curve with Origin.
        FILEOUT2 = [FILE(1:i) 'grat_sf_curv'];
        fileid = [PATHOUT FILEOUT2];
        proffid = fopen(fileid, 'w');
        fprintf(proffid,'SFIn\tFit\tSF\tAvgResp\tStdErr\tSF2\tSpon\n');
        for i=1:length(x_interp)
            fprintf(proffid,'%6.4f\t%6.4f\t', x_interp(i), y_interp(i));
            if (i <= length(px))
                fprintf(proffid,'%6.4f\t%6.2f\t%6.3f\t', px(i), py(i), perr(i));
            else
                fprintf(proffid,'\t\t\t');
            end
            if (i <= 2)
                fprintf(proffid,'%6.4f\t%6.2f\n',null_x(i),null_y(i));
            else
                fprintf(proffid,'\t\n');
            end
        end
        fclose(proffid);
        
        
        %------------------------------------------------------------------------
        %write out all relevant parameters to a cumulative text file, GCD 8/08/01
        %write out one line for each stimu_type for each neuron.
        outfile = [BASE_PATH 'ProtocolSpecific\GratingSpatialFreq\GratingSF_Summary.dat'];
        printflag = 0;
        if (exist(outfile, 'file') == 0)    %file does not yet exist
            printflag = 1;
        end
        fid = fopen(outfile, 'a');
        if (printflag)
            fprintf(fid, 'FILE\t\t\t PrfOR\t PrfSF\t PrfTF\t RFXctr\t RFYctr\t RFDiam\t GrType\t MaxLoc\t MaxRsp\t MinLoc\t MinRsp\t Spont\t SFMI\t\t AnovaP\t SFDI  \t VarTerm\t MaxX \t MinX \t SFHWHM\t SFHW25\t SFHW20\t SFHW15\t SFHW10\t RsqMn\t PMean\t RsqRaw\t Praw\t\t Chi2\t\t ChiP\t\t Brat\t\t Ampl\t\t Peak\t\t SD\t\t Offst\t\t');
            fprintf(fid, '\r\n');
            printflag = 0;
        end
        for j = 1:length(unique_stim_type)
            buff = sprintf('%s\t %5.1f\t %5.2f\t %5.3f\t %5.2f\t %5.2f\t %5.2f\t %5.1f\t %5.4f\t %5.2f\t %5.4f\t %5.2f\t %5.2f\t %5.3f\t %10.8f\t %5.3f\t %6.3f\t %5.2f\t %5.2f\t %7.4f\t %7.4f\t %7.4f\t %7.4f\t %7.4f\t %6.4f\t %10.8f\t %6.4f\t %10.8f\t %6.4f\t %10.8f\t %5.2f\t %5.2f\t %6.4f\t %6.4f\t %6.4f\t', ...
                FILE, data.neuron_params(PREFERRED_ORIENTATION, 1), data.neuron_params(PREFERRED_SPATIAL_FREQ, 1), data.neuron_params(PREFERRED_TEMPORAL_FREQ, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
                unique_stim_type(j), pmax{j}.x, pmax{j}.y, pmin{j}.x, pmin{j}.y, null_rate, SFMI(j), p_value(j), SFDI(j), var_term(j), max(unique_sf), min(unique_sf), sf_width_halfmax(j), sf_width_25_percent(j), sf_width_20_percent(j), sf_width_15_percent(j), sf_width_10_percent(j), stats1(j,1), stats1(j,3), stats2(j,1), stats2(j,3), chi2(j), chiP(j), base_rate(j), amplitude(j), peak(j), st_dev(j), log_offset(j) );
            fprintf(fid, '%s', buff);
            fprintf(fid, '\r\n');
        end
        fclose(fid);
        %------------------------------------------------------------------------
        
    end
    
    %This output code will be used to determine what parameters to use that
    %best describe the differences in the widths of the curves
    output = 0;
    if (output == 1)
        
        %------------------------------------------------------------------------
        %write out all relevant parameters to a cumulative text file, GCD 8/08/01
        %write out one line for each stimu_type for each neuron.
        outfile = [BASE_PATH 'ProtocolSpecific\GratingSpatialFreq\Percentage_Summary.dat'];
        printflag = 0;
        if (exist(outfile, 'file') == 0)    %file does not yet exist
            printflag = 1;
        end
        fid = fopen(outfile, 'a');
        if (printflag)
            fprintf(fid, 'FILE\t\t\t 25\t\t 20\t\t 15\t\t 10\t');
            fprintf(fid, '\r\n');
            printflag = 0;
        end
        for j = 1:length(unique_stim_type)
            buff = sprintf('%s\t %5.5f\t %5.5f\t %5.5f\t %5.5f\t', ...
                FILE, sf_width_25_percent(j), sf_width_20_percent(j), sf_width_15_percent(j), sf_width_10_percent(j));
            fprintf(fid, '%s', buff);
            fprintf(fid, '\r\n');
        end
        fclose(fid);
        %------------------------------------------------------------------------
        
    end  %if (output == 1)
    
%%
%% DHK - 3 Laser Conditions -> two experiments
%%
elseif length(intersect(unique(data.laser_values),[0 1]))==2
    sf = data.gratings_params(GRAT_SPATIAL_FREQ,:,PATCH1);
    null_trials = logical( (sf == data.one_time_params(NULL_VALUE)) );
    stim_type = data.gratings_params(GRAT_TYPE,:,PATCH1);
    unique_stim_type = munique(stim_type(~null_trials)');
    spike_rates = data.spike_rates(SpikeChan, :);
    unique_sf = munique(sf(~null_trials)')
    trials = 1:length(sf);		% a vector of trial indices
    select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );
    
    %% Stimuli Only
    ind_laser0=find(data.laser_values==0); % Stimuli only
    figure(2);
    set(gcf,'PaperPosition', [.51 .51 20 27.5], 'Position', [450 50 500 573], 'Name', 'Grating Spatial Frequency Tuning Curve');
    subplot(3, 1, 3);
    
    for j = 1:length(unique_stim_type);
        
        grating_select = logical( (stim_type == unique_stim_type(j)) );
        
%         plot_x=sf(~null_trials & select_trials & grating_select);
%         plot_y=spike_rates(~null_trials & select_trials & grating_select);
        plot_x=sf(ind_laser0);
        plot_y=spike_rates(ind_laser0);
    
        [px, py, perr] = PlotTuningCurve(plot_x', plot_y', symbols{j}, lines{j}, 1, 0);
        unique_sfreqs=px; %store category groups for ANOVAs
        hold on;
        plot(plot_x, plot_y, 'k.');
        errorbar(px, py, perr, perr, symbols{j});
        hold on;
        
        [SFDI(j), var_term(j)] = Compute_DDI(plot_x, plot_y);
        
        means = [px py];
        raw = [plot_x' plot_y'];
        [pars] = SF_loggaussfit(means, raw); %last arg=0 => only allow positive going fit
        x_interp = (px(1)):0.001:(px(length(px)));
        y_interp = SF_loggaussfunc(x_interp, pars);
        plot(x_interp, y_interp, 'k-');
        null_x = [min(px) max(px)];
        null_resp = data.spike_rates(SpikeChan, null_trials & select_trials);
        null_rate = mean(null_resp);
        null_y = [null_rate null_rate];
        hold on;
        plot(null_x, null_y, 'k--');
        hold off;
        
        yl=ylim;
        ylim([0 yl(2)]); %set lower limit of the Y axis to zero
        xlabel('Grating Spatial Frequency (cyc/deg)');
        ylabel('Response (spikes/sec)');
        y_fit = SF_loggaussfunc(px, pars);
        y_fit(y_fit<0) = 0;
        y_fit = [ones(length(y_fit), 1) y_fit];
        [b, bint, r, rint, stats1(j,:)] = regress(py, y_fit);
        
        y_fit_raw = SF_loggaussfunc(plot_x', pars);
        y_fit_raw(y_fit_raw<0) = 0;
        y_fit_raw = [ones(length(y_fit_raw), 1) y_fit_raw];
        [b, bint, r, rint, stats2(j,:)] = regress(plot_y', y_fit_raw);
        p_value(j) = spk_anova(plot_y, plot_x, unique_sfreqs);
        base_rate(j) = pars(1);
        amplitude(j) = pars(2);
        peak(j) = pars(3);
        st_dev(j) = pars(4);
        log_offset(j) = pars(5);
        max_rate(j) = base_rate(j) + amplitude(j);
        SFMI(j) = 1 - (base_rate - null_rate)/(max_rate - null_rate);
        [chi2(j), chiP(j)] = Chi2_Test(plot_x, plot_y, 'SF_loggaussfunc', pars, length(pars));
        [pmax{j}.y max_i] = max(y_interp);
        pmax{j}.x = x_interp(max_i);
        [pmin{j}.y min_i] = min(y_interp);
        pmin{j}.x = x_interp(min_i);
        
        index = max_i;
        ysearch = pmax{j}.y;
        while (ysearch > (pmin{j}.y + (pmax{j}.y-pmin{j}.y)/2))
            ysearch = y_interp(index);
            if (index == length(y_interp))
                break;
            else
                index = index + 1;
            end
        end
        sf_width_halfmax(j) = x_interp(index) - pmax{j}.x;
        index = max_i;
        ysearch = pmax{j}.y;
        while (ysearch > (pmin{j}.y + (pmax{j}.y-pmin{j}.y)*0.25))
            ysearch = y_interp(index);
            if (index == length(y_interp))
                sf_width_25_percent(j) = 9999;
            else
                index = index + 1;
            end
        end
        sf_width_25_percent(j) = x_interp(index) - pmax{j}.x;
        
        index = max_i;
        ysearch = pmax{j}.y;
        while (ysearch > (pmin{j}.y + (pmax{j}.y-pmin{j}.y)*0.20))
            ysearch = y_interp(index);
            if (index == length(y_interp))
                sf_width_20_percent(j) = 9999;
            else
                index = index + 1;
            end
        end
        sf_width_20_percent(j) = x_interp(index) - pmax{j}.x;
        
        index = max_i;
        ysearch = pmax{j}.y;
        while (ysearch > (pmin{j}.y + (pmax{j}.y-pmin{j}.y)*0.15))
            ysearch = y_interp(index);
            if (index == length(y_interp))
                sf_width_15_percent(j) = 9999;
            else
                index = index + 1;
            end
        end
        sf_width_15_percent(j) = x_interp(index) - pmax{j}.x;
        
        index = max_i;
        ysearch = pmax{j}.y;
        while (ysearch > (pmin{j}.y + (pmax{j}.y-pmin{j}.y)*0.10))
            ysearch = y_interp(index);
            if (index == length(y_interp))
                sf_width_10_percent(j) = 9999;
            else
                index = index + 1;
            end
        end
        sf_width_10_percent(j) = x_interp(index) - pmax{j}.x;
    end
    
    %now, print out some useful information in the upper subplot
    subplot(3, 1, 1);
    PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    PrintSpatialFreqData(p_value, base_rate, amplitude, peak, st_dev, log_offset, max_rate, stats1, stats2, chi2, chiP, SFDI, sf_width_halfmax);
    
    output = 1;
    if (output == 1)        
        i = size(PATH,2) - 1;
        while PATH(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
            i = i - 1;
        end
        PATHOUT = [PATH(1:i) 'Analysis\Tuning\'];
        i = size(FILE,2) - 1;
        while FILE(i) ~='.'
            i = i - 1;
        end
        FILEOUT2 = [FILE(1:i) 'grat_sf_curv'];
        fileid = [PATHOUT FILEOUT2];
        proffid = fopen(fileid, 'w');
        fprintf(proffid,'SFIn\tFit\tSF\tAvgResp\tStdErr\tSF2\tSpon\n');
        for i=1:length(x_interp)
            fprintf(proffid,'%6.4f\t%6.4f\t', x_interp(i), y_interp(i));
            if (i <= length(px))
                fprintf(proffid,'%6.4f\t%6.2f\t%6.3f\t', px(i), py(i), perr(i));
            else
                fprintf(proffid,'\t\t\t');
            end
            if (i <= 2)
                fprintf(proffid,'%6.4f\t%6.2f\n',null_x(i),null_y(i));
            else
                fprintf(proffid,'\t\n');
            end
        end
        fclose(proffid);
        outfile = [BASE_PATH 'ProtocolSpecific\GratingSpatialFreq\GratingSF_Summary.dat'];
        printflag = 0;
        if (exist(outfile, 'file') == 0)    %file does not yet exist
            printflag = 1;
        end
        fid = fopen(outfile, 'a');
        if (printflag)
            fprintf(fid, 'FILE\t\t\t PrfOR\t PrfSF\t PrfTF\t RFXctr\t RFYctr\t RFDiam\t GrType\t MaxLoc\t MaxRsp\t MinLoc\t MinRsp\t Spont\t SFMI\t\t AnovaP\t SFDI  \t VarTerm\t MaxX \t MinX \t SFHWHM\t SFHW25\t SFHW20\t SFHW15\t SFHW10\t RsqMn\t PMean\t RsqRaw\t Praw\t\t Chi2\t\t ChiP\t\t Brat\t\t Ampl\t\t Peak\t\t SD\t\t Offst\t\t');
            fprintf(fid, '\r\n');
            printflag = 0;
        end
        FILE2=strcat(FILE,'ViS');
        for j = 1:length(unique_stim_type)
            buff = sprintf('%s\t %5.1f\t %5.2f\t %5.3f\t %5.2f\t %5.2f\t %5.2f\t %5.1f\t %5.4f\t %5.2f\t %5.4f\t %5.2f\t %5.2f\t %5.3f\t %10.8f\t %5.3f\t %6.3f\t %5.2f\t %5.2f\t %7.4f\t %7.4f\t %7.4f\t %7.4f\t %7.4f\t %6.4f\t %10.8f\t %6.4f\t %10.8f\t %6.4f\t %10.8f\t %5.2f\t %5.2f\t %6.4f\t %6.4f\t %6.4f\t', ...
                FILE2, data.neuron_params(PREFERRED_ORIENTATION, 1), data.neuron_params(PREFERRED_SPATIAL_FREQ, 1), data.neuron_params(PREFERRED_TEMPORAL_FREQ, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
                unique_stim_type(j), pmax{j}.x, pmax{j}.y, pmin{j}.x, pmin{j}.y, null_rate, SFMI(j), p_value(j), SFDI(j), var_term(j), max(unique_sf), min(unique_sf), sf_width_halfmax(j), sf_width_25_percent(j), sf_width_20_percent(j), sf_width_15_percent(j), sf_width_10_percent(j), stats1(j,1), stats1(j,3), stats2(j,1), stats2(j,3), chi2(j), chiP(j), base_rate(j), amplitude(j), peak(j), st_dev(j), log_offset(j) );
            fprintf(fid, '%s', buff);
            fprintf(fid, '\r\n');
        end
        fclose(fid);
    end
    output = 0;
    if (output == 1)outfile = [BASE_PATH 'ProtocolSpecific\GratingSpatialFreq\Percentage_Summary.dat'];
        printflag = 0;
        if (exist(outfile, 'file') == 0)    %file does not yet exist
            printflag = 1;
        end
        fid = fopen(outfile, 'a');
        if (printflag)
            fprintf(fid, 'FILE\t\t\t 25\t\t 20\t\t 15\t\t 10\t');
            fprintf(fid, '\r\n');
            printflag = 0;
        end
        for j = 1:length(unique_stim_type)
            FILE2=strcat(FILE,'ViS');
            buff = sprintf('%s\t %5.5f\t %5.5f\t %5.5f\t %5.5f\t', ...
                FILE2, sf_width_25_percent(j), sf_width_20_percent(j), sf_width_15_percent(j), sf_width_10_percent(j));
            fprintf(fid, '%s', buff);
            fprintf(fid, '\r\n');
        end
        fclose(fid);
    end
    
    %% Stimuli + Laser
    ind_laser1=find(data.laser_values==1); % Stimuli + Laser
    subplot(3, 1, 3);
    
    for j = 1:length(unique_stim_type);
        
        grating_select = logical( (stim_type == unique_stim_type(j)) );
        
%         plot_x=sf(~null_trials & select_trials & grating_select);
%         plot_y=spike_rates(~null_trials & select_trials & grating_select);
        plot_x=sf(ind_laser1);
        plot_y=spike_rates(ind_laser1);
    
        [px, py, perr] = PlotTuningCurve(plot_x', plot_y', symbols{j}, lines{j}, 1, 0);
        unique_sfreqs=px; %store category groups for ANOVAs
        hold on;
        plot(plot_x, plot_y, 'r.');
        errorbar(px, py, perr, perr, symbols{j});
        hold on;
        
        [SFDI(j), var_term(j)] = Compute_DDI(plot_x, plot_y);
        
        means = [px py];
        raw = [plot_x' plot_y'];
        [pars] = SF_loggaussfit(means, raw); %last arg=0 => only allow positive going fit
        x_interp = (px(1)):0.001:(px(length(px)));
        y_interp = SF_loggaussfunc(x_interp, pars);
        plot(x_interp, y_interp, 'r-');
        null_x = [min(px) max(px)];
        null_resp = data.spike_rates(SpikeChan, null_trials & select_trials);
        null_rate = mean(null_resp);
        null_y = [null_rate null_rate];
        hold on;
        plot(null_x, null_y, 'r--');
        hold off;
        
        yl=ylim;
        ylim([0 yl(2)]); %set lower limit of the Y axis to zero
        xlabel('Grating Spatial Frequency (cyc/deg)');
        ylabel('Response (spikes/sec)');
        y_fit = SF_loggaussfunc(px, pars);
        y_fit(y_fit<0) = 0;
        y_fit = [ones(length(y_fit), 1) y_fit];
        [b, bint, r, rint, stats1(j,:)] = regress(py, y_fit);
        
        y_fit_raw = SF_loggaussfunc(plot_x', pars);
        y_fit_raw(y_fit_raw<0) = 0;
        y_fit_raw = [ones(length(y_fit_raw), 1) y_fit_raw];
        [b, bint, r, rint, stats2(j,:)] = regress(plot_y', y_fit_raw);
        p_value(j) = spk_anova(plot_y, plot_x, unique_sfreqs);
        base_rate(j) = pars(1);
        amplitude(j) = pars(2);
        peak(j) = pars(3);
        st_dev(j) = pars(4);
        log_offset(j) = pars(5);
        max_rate(j) = base_rate(j) + amplitude(j);
        SFMI(j) = 1 - (base_rate - null_rate)/(max_rate - null_rate);
        [chi2(j), chiP(j)] = Chi2_Test(plot_x, plot_y, 'SF_loggaussfunc', pars, length(pars));
        [pmax{j}.y max_i] = max(y_interp);
        pmax{j}.x = x_interp(max_i);
        [pmin{j}.y min_i] = min(y_interp);
        pmin{j}.x = x_interp(min_i);
        
        index = max_i;
        ysearch = pmax{j}.y;
        while (ysearch > (pmin{j}.y + (pmax{j}.y-pmin{j}.y)/2))
            ysearch = y_interp(index);
            if (index == length(y_interp))
                break;
            else
                index = index + 1;
            end
        end
        sf_width_halfmax(j) = x_interp(index) - pmax{j}.x;
        index = max_i;
        ysearch = pmax{j}.y;
        while (ysearch > (pmin{j}.y + (pmax{j}.y-pmin{j}.y)*0.25))
            ysearch = y_interp(index);
            if (index == length(y_interp))
                sf_width_25_percent(j) = 9999;
            else
                index = index + 1;
            end
        end
        sf_width_25_percent(j) = x_interp(index) - pmax{j}.x;
        
        index = max_i;
        ysearch = pmax{j}.y;
        while (ysearch > (pmin{j}.y + (pmax{j}.y-pmin{j}.y)*0.20))
            ysearch = y_interp(index);
            if (index == length(y_interp))
                sf_width_20_percent(j) = 9999;
            else
                index = index + 1;
            end
        end
        sf_width_20_percent(j) = x_interp(index) - pmax{j}.x;
        
        index = max_i;
        ysearch = pmax{j}.y;
        while (ysearch > (pmin{j}.y + (pmax{j}.y-pmin{j}.y)*0.15))
            ysearch = y_interp(index);
            if (index == length(y_interp))
                sf_width_15_percent(j) = 9999;
            else
                index = index + 1;
            end
        end
        sf_width_15_percent(j) = x_interp(index) - pmax{j}.x;
        
        index = max_i;
        ysearch = pmax{j}.y;
        while (ysearch > (pmin{j}.y + (pmax{j}.y-pmin{j}.y)*0.10))
            ysearch = y_interp(index);
            if (index == length(y_interp))
                sf_width_10_percent(j) = 9999;
            else
                index = index + 1;
            end
        end
        sf_width_10_percent(j) = x_interp(index) - pmax{j}.x;
    end
    
    %now, print out some useful information in the upper subplot
    subplot(3, 1, 2);
   PrintSpatialFreqData2(p_value, base_rate, amplitude, peak, st_dev, log_offset, max_rate, stats1, stats2, chi2, chiP, SFDI, sf_width_halfmax);
    
    output = 1;
    if (output == 1)        
        i = size(PATH,2) - 1;
        while PATH(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
            i = i - 1;
        end
        PATHOUT = [PATH(1:i) 'Analysis\Tuning\'];
        i = size(FILE,2) - 1;
        while FILE(i) ~='.'
            i = i - 1;
        end
        FILEOUT2 = [FILE(1:i) 'grat_sf_curv'];
        fileid = [PATHOUT FILEOUT2];
        proffid = fopen(fileid, 'w');
        fprintf(proffid,'SFIn\tFit\tSF\tAvgResp\tStdErr\tSF2\tSpon\n');
        for i=1:length(x_interp)
            fprintf(proffid,'%6.4f\t%6.4f\t', x_interp(i), y_interp(i));
            if (i <= length(px))
                fprintf(proffid,'%6.4f\t%6.2f\t%6.3f\t', px(i), py(i), perr(i));
            else
                fprintf(proffid,'\t\t\t');
            end
            if (i <= 2)
                fprintf(proffid,'%6.4f\t%6.2f\n',null_x(i),null_y(i));
            else
                fprintf(proffid,'\t\n');
            end
        end
        fclose(proffid);
        outfile = [BASE_PATH 'ProtocolSpecific\GratingSpatialFreq\GratingSF_Summary.dat'];
        printflag = 0;
        if (exist(outfile, 'file') == 0)    %file does not yet exist
            printflag = 1;
        end
        fid = fopen(outfile, 'a');
        if (printflag)
            fprintf(fid, 'FILE\t\t\t PrfOR\t PrfSF\t PrfTF\t RFXctr\t RFYctr\t RFDiam\t GrType\t MaxLoc\t MaxRsp\t MinLoc\t MinRsp\t Spont\t SFMI\t\t AnovaP\t SFDI  \t VarTerm\t MaxX \t MinX \t SFHWHM\t SFHW25\t SFHW20\t SFHW15\t SFHW10\t RsqMn\t PMean\t RsqRaw\t Praw\t\t Chi2\t\t ChiP\t\t Brat\t\t Ampl\t\t Peak\t\t SD\t\t Offst\t\t');
            fprintf(fid, '\r\n');
            printflag = 0;
        end
            FILE3=strcat(FILE,'V+L');
        for j = 1:length(unique_stim_type)
            buff = sprintf('%s\t %5.1f\t %5.2f\t %5.3f\t %5.2f\t %5.2f\t %5.2f\t %5.1f\t %5.4f\t %5.2f\t %5.4f\t %5.2f\t %5.2f\t %5.3f\t %10.8f\t %5.3f\t %6.3f\t %5.2f\t %5.2f\t %7.4f\t %7.4f\t %7.4f\t %7.4f\t %7.4f\t %6.4f\t %10.8f\t %6.4f\t %10.8f\t %6.4f\t %10.8f\t %5.2f\t %5.2f\t %6.4f\t %6.4f\t %6.4f\t', ...
                FILE3, data.neuron_params(PREFERRED_ORIENTATION, 1), data.neuron_params(PREFERRED_SPATIAL_FREQ, 1), data.neuron_params(PREFERRED_TEMPORAL_FREQ, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
                unique_stim_type(j), pmax{j}.x, pmax{j}.y, pmin{j}.x, pmin{j}.y, null_rate, SFMI(j), p_value(j), SFDI(j), var_term(j), max(unique_sf), min(unique_sf), sf_width_halfmax(j), sf_width_25_percent(j), sf_width_20_percent(j), sf_width_15_percent(j), sf_width_10_percent(j), stats1(j,1), stats1(j,3), stats2(j,1), stats2(j,3), chi2(j), chiP(j), base_rate(j), amplitude(j), peak(j), st_dev(j), log_offset(j) );
            fprintf(fid, '%s', buff);
            fprintf(fid, '\r\n');
        end
        fclose(fid);
    end
    output = 0;
    if (output == 1)outfile = [BASE_PATH 'ProtocolSpecific\GratingSpatialFreq\Percentage_Summary.dat'];
        printflag = 0;
        if (exist(outfile, 'file') == 0)    %file does not yet exist
            printflag = 1;
        end
        fid = fopen(outfile, 'a');
        if (printflag)
            fprintf(fid, 'FILE\t\t\t 25\t\t 20\t\t 15\t\t 10\t');
            fprintf(fid, '\r\n');
            printflag = 0;
        end
        for j = 1:length(unique_stim_type)
            FILE3=strcat(FILE,'V+L');
            buff = sprintf('%s\t %5.5f\t %5.5f\t %5.5f\t %5.5f\t', ...
                FILE3, sf_width_25_percent(j), sf_width_20_percent(j), sf_width_15_percent(j), sf_width_10_percent(j));
            fprintf(fid, '%s', buff);
            fprintf(fid, '\r\n');
        end
        fclose(fid);
    end
%%
end
return;