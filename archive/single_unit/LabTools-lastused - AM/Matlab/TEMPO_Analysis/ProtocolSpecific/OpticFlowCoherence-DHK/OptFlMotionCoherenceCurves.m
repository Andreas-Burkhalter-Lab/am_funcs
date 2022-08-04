%-----------------------------------------------------------------------------------------------------------------------
%-- MotionCoherenceCurves.m -- Plots response vs motion coherence for preferred
%and null directions of motion, GCD, 3/22/06
%-----------------------------------------------------------------------------------------------------------------------

%%%% 4/27/16 AM added conditional when writing outfile to write into the correct
%%%% directory if on a specific computer; functionality should be
%%%% unchanged if using on any lab computer
%%%% 4/27/16 AM edited formatting of output file columns 
% 3/14/16 AM corrected typos, updated file paths
%%%% 8/23/16 AM added optional 'skip_trials' argin - these trials will not be analyzed;
%%%%    trials not included within BegTrial-EndTrial will still be skipped

function OptFlMotionCoherenceCurves(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, skip_trials)

TEMPO_Defs;
Path_Defs;

symbols = {'ko' 'k*' 'go' 'mo' 'b*' 'r*' 'g*' 'c*'};
lines = {'k-' 'k--' 'g-' 'm-' 'b--' 'r--' 'g--' 'c--'};

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

%%%% AM added 8/29/16
if exist('skip_trials','var')
    select_trials(skip_trials) = false;
end


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
    plot(plot_x, plot_y, 'k.');
    errorbar(px, py, perr, perr, symbols{i});
    hold on;

    %Fit the data with a straight line and plot the line
    p(i,:) = polyfit(px, py, 1);
    x_interp = [min(px):0.11:ceil(max(px))];
    y_interp = x_interp.*p(i,1) + p(i,2);
    plot(x_interp, y_interp, lines{i})
    
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
plot(null_x, null_y, 'k--');
hold off;

yl = ylim;
ylim([0 yl(2)]);	% set the lower limit of the Y axis to zero
xlabel('coherence');
ylabel('Response (spikes/sec)');


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
    
    %%%% 4/27/16 AM added conditional below to write into the correct
    %%%% directory if on a specific computer; if on any other lab computer
    %%%% such as 'BURK-LV', the outfile will be written to the usual directory
    if strcmp(getenv('COMPUTERNAME'),'ANDREW') %%% Andrew's MSI laptop name
        outfile = 'C:\Users\AM\Documents\Lab\recordings\OpticFlowCoherence_Summary.dat';
    else
        outfile = [BASE_PATH 'ProtocolSpecific\OpticFlowCoherence-DHK\OpticFlowCoherence_Summary.dat'];
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
        %%% original line: fprintf(fid, 'FILE\t\t\t PrfOR\t PrfSF\t PrfTF\t RFXctr\t RFYctr\t RFDiam\t Direction\t\t Slope\t Int\t\t RsqMn\t PMn\t\t RsqRaw\t PRaw\t\t ');
        fprintf(fid, 'FILE\t PrfOR\t PrfSF\t PrfTF\t RFXctr\t RFYctr\t RFDiam\t Direction\t Slope\t Int\t RsqMn\t PMn\t RsqRaw\t PRaw\t ');
        fprintf(fid, '\r\n');
        printflag = 0;
    end
    for j = 1:length(unique_direc)
        %%%% 4/27/16 AM edited following line - removed all repeated
        %%%% instances of '\t', which added extra unnecessary columns to
        %%%% output file
        %%% original line: buff = sprintf('%s\t %5.1f\t %5.2f\t %5.3f\t %5.2f\t %5.2f\t %5.2f\t %4.4f\t\t %5.4f\t %5.4f\t %5.4f\t %5.4f\t %5.4f\t %5.4f\t ', ...
        %%%                    FILE, data.neuron_params(PREFERRED_ORIENTATION, 1), data.neuron_params(PREFERRED_SPATIAL_FREQ, 1), data.neuron_params(PREFERRED_TEMPORAL_FREQ, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
        %%%                    unique_direc(j), p(j,1), p(j,2), Rsq_means(j), P_means(j), Rsq_raw(j), P_raw(j) );
        buff = sprintf('%s\t %5.1f\t %5.2f\t %5.3f\t %5.2f\t %5.2f\t %5.2f\t %4.4f\t %5.4f\t %5.4f\t %5.4f\t %5.4f\t %5.4f\t %5.4f\t ', ...
            FILE, data.neuron_params(PREFERRED_ORIENTATION, 1), data.neuron_params(PREFERRED_SPATIAL_FREQ, 1), data.neuron_params(PREFERRED_TEMPORAL_FREQ, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
            unique_direc(j), p(j,1), p(j,2), Rsq_means(j), P_means(j), Rsq_raw(j), P_raw(j) );
        fprintf(fid, '%s', buff);
        fprintf(fid, '\r\n');
    end
    fclose(fid);
    %------------------------------------------------------------------------

end



return;