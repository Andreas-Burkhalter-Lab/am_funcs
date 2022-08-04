%-----------------------------------------------------------------------------------------------------------------------
%-- MotionCoherenceCurves.m -- Plots response vs motion coherence for preferred
%and null directions of motion, GCD, 3/22/06
%-----------------------------------------------------------------------------------------------------------------------
function MotionCoherenceCurves(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;
Path_Defs;

symbols = {'ko' 'k*' 'go' 'mo' 'b*' 'r*' 'g*' 'c*'};
lines = {'k-' 'k--' 'g-' 'm-' 'b--' 'r--' 'g--' 'c--'};

if length(intersect(unique(data.laser_values),[0 1]))==1
    %% DHK - 2 Laser Conditon -> one experiments
    %get the column of values of directions in the dots_params matrix
    coher = data.dots_params(DOTS_COHER,:,PATCH1);
    direc = data.dots_params(DOTS_DIREC,:,PATCH1);
    
    %now, get the firing rates for all the trials
    spike_rates = data.spike_rates(SpikeChan, :);
    
    %get indices of any NULL conditions (for measuring spontaneous activity
    null_trials = logical( (coher == data.one_time_params(NULL_VALUE)) );
    
    %now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
    trials = 1:length(coher);		% a vector of trial indices
    select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );
    
    unique_coher = munique(coher(~null_trials)');
    unique_direc = munique(direc(~null_trials)');
    
    figure;
    set(gcf,'PaperPosition', [.51 .51 20 27.5], 'Position', [450 50 500 573], 'Name', 'Coherence Curves');
    subplot(2, 1, 2);
    
    for i=1:length(unique_direc)	%for each different motion direction (pref and null), plot a separate tuning curve
        direc_select = logical( (direc == unique_direc(i)) );
        
        plot_x = coher(direc_select & ~null_trials & select_trials);
        plot_y = spike_rates(direc_select & ~null_trials & select_trials);
        
        hold on;
        %NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
        [px, py, perr] = PlotTuningCurve(plot_x', plot_y', symbols{i}, 0, 0, 0);
        
        %Plot the data
        hold on;
        plot(log(plot_x), plot_y, 'k.');
        errorbar(log(px), py, perr, perr, symbols{i});
        hold on;
        
        %Fit the data with a straight line and plot the line
        p(i,:) = polyfit(px, py, 1);
        x_interp = [min(px):0.11:ceil(max(px))];
        y_interp = x_interp.*p(i,1) + p(i,2);
        plot(log(x_interp), y_interp, lines{i})
        
        %Compute R^2 for both means and raw values
        y_fit = px.*p(i,1) + p(i,2);
        y_fit(y_fit<0) = 0;
        %add a column of ones to yfit to make regress happy
        y_fit = [ones(length(y_fit), 1) y_fit];
        [b, bint, r, rint, stats1(i,:)] = regress(py, y_fit);
        
        y_fit_raw = plot_x'.*p(i,1) + p(i,2);
        y_fit_raw(y_fit_raw<0) = 0;
        y_fit_raw = [ones(length(y_fit_raw), 1) y_fit_raw];
        [b, bint, r, rint, stats2(i,:)] = regress(plot_y', y_fit_raw);
        
        
        
    end
    
    Rsq_means = stats1(:,1);
    P_means = stats1(:,3);
    Rsq_raw = stats2(:,1);
    P_raw = stats2(:,3);
    
    %now, get the firing rate for NULL condition trials and add spontaneous rate to plot
    null_x = [min(px) max(px)];
    null_rate = mean(data.spike_rates(SpikeChan, null_trials & select_trials));
    null_y = [null_rate null_rate];
    hold on;
    plot(log(null_x), null_y, 'k--');
    hold off;
    
    yl = ylim; %Ji change YLim to ylim10/5/2015
    ylim([0 yl(2)]);	% set the lower limit of the Y axis to zero
    xlabel('log(percent coherence)');
    ylabel('Response (spikes/sec)');%JI change XL, YL to xl, yl
    
    
    %now, print out some useful information in the upper subplot
    subplot(2, 1, 1);
    PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    
    PrintCoherenceData(unique_direc, p, Rsq_means, P_means, Rsq_raw, P_raw);
    %(p_value, base_rate, amplitude, peak, st_dev, log_offset, max_rate, stats1, stats2, chi2, chiP, SFDI, sf_width_halfmax);
    
    output = 1;
    if (output == 1)
        
        %------------------------------------------------------------------------
        %write out all relevant parameters to a cumulative text file, GCD 8/08/01
        %write out one line for each stimu_type for each neuron.
        outfile = [BASE_PATH 'ProtocolSpecific\Coherence\Coherence_Summary.dat'];
        printflag = 0;
        if (exist(outfile, 'file') == 0)    %file does not yet exist
            printflag = 1;
        end
        fid = fopen(outfile, 'a');
        if (printflag)
            fprintf(fid, 'FILE\t\t\t PrfOR\t PrfSF\t PrfTF\t RFXctr\t RFYctr\t RFDiam\t Direction\t\t Slope\t Int\t\t RsqMn\t PMn\t\t RsqRaw\t PRaw\t\t ');
            fprintf(fid, '\r\n');
            printflag = 0;
        end
        for j = 1:length(unique_direc)
            buff = sprintf('%s\t %5.1f\t %5.2f\t %5.3f\t %5.2f\t %5.2f\t %5.2f\t %4.4f\t\t %5.4f\t %5.4f\t %5.4f\t %5.4f\t %5.4f\t %5.4f\t ', ...
                FILE, data.neuron_params(PREFERRED_ORIENTATION, 1), data.neuron_params(PREFERRED_SPATIAL_FREQ, 1), data.neuron_params(PREFERRED_TEMPORAL_FREQ, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
                unique_direc(j), p(j,1), p(j,2), Rsq_means(j), P_means(j), Rsq_raw(j), P_raw(j) );
            fprintf(fid, '%s', buff);
            fprintf(fid, '\r\n');
        end
        fclose(fid);
        %------------------------------------------------------------------------
        
    end
    
%%
%% DHK - 3 Laser Conditions -> two experiments
%%
elseif length(intersect(unique(data.laser_values),[0 1]))==2
    %get the column of values of directions in the dots_params matrix
    coher = data.dots_params(DOTS_COHER,:,PATCH1);
    direc = data.dots_params(DOTS_DIREC,:,PATCH1);
    
    %now, get the firing rates for all the trials
    spike_rates = data.spike_rates(SpikeChan, :);
    
    %get indices of any NULL conditions (for measuring spontaneous activity
    null_trials = logical( (coher == data.one_time_params(NULL_VALUE)) );
    
    %now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
    trials = 1:length(coher);		% a vector of trial indices
    select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );
    
    unique_coher = munique(coher(~null_trials)');
    unique_direc = munique(direc(~null_trials)');
    
    %% Stimuli Only
    ind_laser0=find(data.laser_values==0); % Stimuli only
    figure(2)
    set(gcf,'PaperPosition', [.51 .51 20 27.5], 'Position', [450 50 500 573], 'Name', 'Coherence Curves');
    subplot(3, 1, 3);
    for i=1:length(unique_direc)	%for each different motion direction (pref and null), plot a separate tuning curve
        direc_select = logical( (direc == unique_direc(i)) );
        
%         plot_x = coher(direc_select & ~null_trials & select_trials);
%         plot_y = spike_rates(direc_select & ~null_trials & select_trials);
         plot_x=coher(ind_laser0);
        plot_y=spike_rates(ind_laser0);
        
        hold on;
        [px, py, perr] = PlotTuningCurve(plot_x', plot_y', symbols{i}, 0, 0, 0);
        hold on;
        plot(log(plot_x), plot_y, 'k.');
        errorbar(log(px), py, perr, perr, symbols{i});
        hold on;
        p(i,:) = polyfit(px, py, 1);
        x_interp = [min(px):0.11:ceil(max(px))];
        y_interp = x_interp.*p(i,1) + p(i,2);
        plot(log(x_interp), y_interp, lines{i})
        y_fit = px.*p(i,1) + p(i,2);
        y_fit(y_fit<0) = 0;
        y_fit = [ones(length(y_fit), 1) y_fit];
        [b, bint, r, rint, stats1(i,:)] = regress(py, y_fit);
        y_fit_raw = plot_x'.*p(i,1) + p(i,2);
        y_fit_raw(y_fit_raw<0) = 0;
        y_fit_raw = [ones(length(y_fit_raw), 1) y_fit_raw];
        [b, bint, r, rint, stats2(i,:)] = regress(plot_y', y_fit_raw);
    end
    Rsq_means = stats1(:,1);
    P_means = stats1(:,3);
    Rsq_raw = stats2(:,1);
    P_raw = stats2(:,3);
    null_x = [min(px) max(px)];
    null_rate = mean(data.spike_rates(SpikeChan, null_trials & select_trials));
    null_y = [null_rate null_rate];
    hold on;
    plot(log(null_x), null_y, 'k--');
    hold off;
    yl = ylim; %Ji change YLim to ylim10/5/2015
    ylim([0 yl(2)]);	% set the lower limit of the Y axis to zero
    xlabel('log(percent coherence)');
    ylabel('Response (spikes/sec)');%JI change XL, YL to xl, yl
    subplot(3, 1, 1);
    PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
    PrintCoherenceData(unique_direc, p, Rsq_means, P_means, Rsq_raw, P_raw);
    output = 1;
    if (output == 1)
        outfile = [BASE_PATH 'ProtocolSpecific\Coherence\Coherence_Summary.dat'];
        printflag = 0;
        if (exist(outfile, 'file') == 0)    %file does not yet exist
            printflag = 1;
        end
        fid = fopen(outfile, 'a');
        if (printflag)
            fprintf(fid, 'FILE\t\t\t PrfOR\t PrfSF\t PrfTF\t RFXctr\t RFYctr\t RFDiam\t Direction\t\t Slope\t Int\t\t RsqMn\t PMn\t\t RsqRaw\t PRaw\t\t ');
            fprintf(fid, '\r\n');
            printflag = 0;
        end
        FILE2=strcat(FILE,'ViS');
        for j = 1:length(unique_direc)
            buff = sprintf('%s\t %5.1f\t %5.2f\t %5.3f\t %5.2f\t %5.2f\t %5.2f\t %4.4f\t\t %5.4f\t %5.4f\t %5.4f\t %5.4f\t %5.4f\t %5.4f\t ', ...
                FILE2, data.neuron_params(PREFERRED_ORIENTATION, 1), data.neuron_params(PREFERRED_SPATIAL_FREQ, 1), data.neuron_params(PREFERRED_TEMPORAL_FREQ, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
                unique_direc(j), p(j,1), p(j,2), Rsq_means(j), P_means(j), Rsq_raw(j), P_raw(j) );
            fprintf(fid, '%s', buff);
            fprintf(fid, '\r\n');
        end
        fclose(fid);
    end
    %% Stimuli + Laser
    ind_laser1=find(data.laser_values==1); % Stimuli + Laser
    
    figure(2)
    subplot(3, 1, 3);
    
    for i=1:length(unique_direc)	%for each different motion direction (pref and null), plot a separate tuning curve
        direc_select = logical( (direc == unique_direc(i)) );
        
%         plot_x = coher(direc_select & ~null_trials & select_trials);
%         plot_y = spike_rates(direc_select & ~null_trials & select_trials);
         plot_x=coher(ind_laser1);
        plot_y=spike_rates(ind_laser1);
        
        hold on;
        [px, py, perr] = PlotTuningCurve(plot_x', plot_y', symbols{i}, 0, 0, 0);
        hold on;
        plot(log(plot_x), plot_y, 'r.');
        errorbar(log(px), py, perr, perr, symbols{i});
        hold on;
        p(i,:) = polyfit(px, py, 1);
        x_interp = [min(px):0.11:ceil(max(px))];
        y_interp = x_interp.*p(i,1) + p(i,2);
        plot(log(x_interp), y_interp,'r-')
        y_fit = px.*p(i,1) + p(i,2);
        y_fit(y_fit<0) = 0;
        y_fit = [ones(length(y_fit), 1) y_fit];
        [b, bint, r, rint, stats1(i,:)] = regress(py, y_fit);
        y_fit_raw = plot_x'.*p(i,1) + p(i,2);
        y_fit_raw(y_fit_raw<0) = 0;
        y_fit_raw = [ones(length(y_fit_raw), 1) y_fit_raw];
        [b, bint, r, rint, stats2(i,:)] = regress(plot_y', y_fit_raw);
    end
    Rsq_means = stats1(:,1);
    P_means = stats1(:,3);
    Rsq_raw = stats2(:,1);
    P_raw = stats2(:,3);
    null_x = [min(px) max(px)];
    null_rate = mean(data.spike_rates(SpikeChan, null_trials & select_trials));
    null_y = [null_rate null_rate];
    hold on;
    plot(log(null_x), null_y, 'r--');
    hold off;
    yl = ylim; %Ji change YLim to ylim10/5/2015
    ylim([0 yl(2)]);	% set the lower limit of the Y axis to zero
    xlabel('log(percent coherence)');
    ylabel('Response (spikes/sec)');%JI change XL, YL to xl, yl
    subplot(3, 1, 2);
    PrintCoherenceData2(unique_direc, p, Rsq_means, P_means, Rsq_raw, P_raw);
    output = 1;
    if (output == 1)
        outfile = [BASE_PATH 'ProtocolSpecific\Coherence\Coherence_Summary.dat'];
        printflag = 0;
        if (exist(outfile, 'file') == 0)    %file does not yet exist
            printflag = 1;
        end
        fid = fopen(outfile, 'a');
        if (printflag)
            fprintf(fid, 'FILE\t\t\t PrfOR\t PrfSF\t PrfTF\t RFXctr\t RFYctr\t RFDiam\t Direction\t\t Slope\t Int\t\t RsqMn\t PMn\t\t RsqRaw\t PRaw\t\t ');
            fprintf(fid, '\r\n');
            printflag = 0;
        end
        FILE3=strcat(FILE,'V+L');
        for j = 1:length(unique_direc)
            buff = sprintf('%s\t %5.1f\t %5.2f\t %5.3f\t %5.2f\t %5.2f\t %5.2f\t %4.4f\t\t %5.4f\t %5.4f\t %5.4f\t %5.4f\t %5.4f\t %5.4f\t ', ...
                FILE3, data.neuron_params(PREFERRED_ORIENTATION, 1), data.neuron_params(PREFERRED_SPATIAL_FREQ, 1), data.neuron_params(PREFERRED_TEMPORAL_FREQ, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
                unique_direc(j), p(j,1), p(j,2), Rsq_means(j), P_means(j), Rsq_raw(j), P_raw(j) );
            fprintf(fid, '%s', buff);
            fprintf(fid, '\r\n');
        end
        fclose(fid);
    end
end
return;