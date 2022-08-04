%-----------------------------------------------------------------------------------------------------------------------
%-- TFTuningCurve.m -- Plots a TF tuning curve with fit with a Gaussian
%
%--	CMA 04/10/06
%-----------------------------------------------------------------------------------------------------------------------
function TFTuningCurve(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;		%needed for defines like IN_T1_WIN_CD
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01
Path_Defs;

symbols = {'ko' 'k*' 'go' 'mo' 'b*' 'r*' 'g*' 'c*'};
lines = {'k-' 'k--' 'g-' 'm-' 'b--' 'r--' 'g--' 'c--'};


%get the column of values of TFs in gratings_params matrix
tf = data.gratings_params(GRAT_TEMPORAL_FREQ,:,PATCH1);

%get indices of any NULL conditions (for measuring spontaneous activity)
null_trials = logical( (tf == data.one_time_params(NULL_VALUE)) );

%get the column of stimulus type values
stim_type = data.gratings_params(GRAT_TYPE,:,PATCH1);
unique_stim_type = munique(stim_type(~null_trials)');

%now, get the firing rates for all the trials 
spike_rates = data.spike_rates(SpikeChan, :);

unique_tf = munique(tf(~null_trials)')

%now, remove trials that do not fall between BegTrial and EndTrial
trials = 1:length(tf);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

plot_x = tf(~null_trials & select_trials);
plot_y = spike_rates(~null_trials & select_trials);

figure(2);
set(gcf,'PaperPosition', [.51 .51 20 27.5], 'Position', [450 50 500 573], 'Name', 'Grating Temporal Frequency Tuning Curve');
subplot(2, 1, 2);

%NOTE: Inputs to PlotTuningCurve.m must be column vectors, not row vectors
%because of munique().
[px, py, perr] = PlotTuningCurve(plot_x', plot_y', 'ko', 'k-', 1, 0); %NOTE: last arg=0 means just get output, no plot

%keep a copy of the original data before shifting below
px_orig = px;
py_orig = py;

unique_temps = px; % store category groups for ANOVAs

%NOTE: In the directional tuning protocol, this section was used to make
%sure that the peak of the tuning curve was always plotted in the center of
%the graph.  If necessary it could be adapted to this protocol, but for the
%time it is being left out.  Chris Andrews 04/07/06
%
%
%--------------------------------------------------------------------------
%now, shift the px, py, and perr vectors such that the peak of the tuning
%curve is in the middle of the axis range
%alternate explanation:  shift these vectors so that the peak response is
%always in the middle of the vector
% ctr_indx = round(length(px)/2 - rem(length(px),2)) + 1;
% [max_val max_indx] = max(py);
% shift = max_indx - ctr_indx;
% if (shift > 0)
%     px = [ px(shift+1 : length(px)) ; px(1 : shift) + 360];
%     px_orig = [px_orig(shift+1 : length(px_orig)) ; px_orig(1 : shift)];
%     py = [py(shift+1 : length(py)) ; py(1 : shift)];
%     perr = [perr(shift+1 : length(perr)) ; perr(1 : shift)];
% end
% 
% if (shift < 0)
%     px = [ px(shift+1 : length(px)) ; px(1 : shift) + 360];
%     px_orig = [px_orig(shift+1 : length(px_orig)) ; px_orig(1 : shift)];
%     py = [py(shift+1 : length(py)) ; py(1 : shift)];
%     perr = [perr(shift+1 : length(perr)) ; perr(1 : shift)];
% end 
%Section continues, see DirectionTuningCurve.m ...
%--------------------------------------------------------------------------

%plot the data
hold on;
plot(plot_x, plot_y, 'k.')
errorbar(px, py, perr, perr, 'ko');
hold on;

%Compute the discrimination index
[TFDI, var_term] = Compute_DDI(plot_x, plot_y);

%Now, fit the data with a Gaussian curve and plot this as well
means = [px py];
raw = [plot_x' plot_y'];
[pars] = gaussfit(means, raw, 0);  %last arg => only allow positive going fit
x_interp = (px(1)): 0.5 : (px(length(px)));
y_interp = gaussfunc(x_interp, pars);
plot(x_interp, y_interp, 'k-');

%now, get the firing rate for NULL condition trials and add spontaneous
%rate to plot
null_x = [min(px) max(px)];
null_resp = data.spike_rates(SpikeChan, null_trials & select_trials);
null_rate = mean(null_resp);
null_y = [null_rate null_rate]
hold on;
plot(null_x, null_y, 'k--');
hold off;

yl = ylim;
ylim([0 yl(2)]); %set the lower limit of the Y axis to zero
xlabel('Grating Temporal Frequency (cyc/sec)');
ylabel('Response (spikes/sec)');

%Compute R^2 of the fit for both means and raw values
y_fit = gaussfunc(px, pars);
y_fit(y_fit < 0) = 0;
%add a column of ones to yfit to make regress happy
y_fit = [ones(length(y_fit), 1) y_fit];
[b, bint, r, rint, stats1] = regress(py, y_fit);

y_fit_raw = gaussfunc(plot_x', pars);
y_fit_raw(y_fit_raw<0) = 0;
y_fit_raw = [ones(length(y_fit_raw), 1) y_fit_raw];
[b, bint, r, rint, stats2] = regress(plot_y', y_fit_raw);



%Now, print out some useful information on the upper subplot
subplot(2,1,1);
PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
%now print out some useful values
%pmax, pmin, py

%    % do the ANOVA on sqrt(firing rate) a la Prince et al.
%    [p_value(i), MSgroup(i), MSerror(i)] = spk_anova_F(sqrt(plot_y), plot_x, px);
%    avg_resp(i) = mean(plot_y);      


%PrintGratingData(p_value, avg_resp, pmax, pmin, px, null_rate, unique_stim_type, PATH, FILE, DDI, DTI);

%Calculate some metrics and stats then print them in the plot
pref_tempfreq = pars(3);
p_value = spk_anova(plot_y, plot_x, unique_temps);
base_rate = pars(1);
amplitude = pars(2);
max_rate = base_rate + amplitude;
width = sqrt(-(log(.5)))*pars(4)*2*sqrt(2);
DSI = 1 - (base_rate - null_rate)/(max_rate - null_rate);

%Do chi-square goodness of fit test
[chi2, chiP] = Chi2_Test(plot_x, plot_y, 'gaussfunc', pars, length(pars));

%calculate modulation index using sqrt raw responses and subtracting
%spontaneous
DMI = Compute_ModIndex(plot_x, plot_y, null_resp);
%PrintDirectionData(p_value, base_rate, null_rate, amplitude, pref_dir, max_rate, width, DSI, stats1, stats2, chi2, chiP);
PrintGratingData(p_value, base_rate, amplitude, pref_tempfreq, max_rate, width, stats1, stats2, DSI, chi2, chiP);

    %------------------------------------------------------------------------
    %write out all relevant parameters to a cumulative text file, GCD 8/08/01
    %write out one line for each stimu_type for each neuron.
    outfile = [BASE_PATH 'ProtocolSpecific\GratingTemporalFreq\GratingTF_Summary.dat'];
    printflag = 0;
    if (exist(outfile, 'file') == 0)    %file does not yet exist
        printflag = 1;
    end
    fid = fopen(outfile, 'a');
    if (printflag)
        fprintf(fid, 'FILE\t\t\t PrfOR\t PrfSF\t PrfTF\t RFXctr\t RFYctr\t RFDiam\t GrType\t MaxLoc\t MaxRsp\t MinLoc\t MinRsp\t AvgRsp\t Spont\t ModIx\t DTI  \t AnovaP\t DDI  \t VarTerm\t MSgrp\t MSerror\t MaxX \t MinX \t');
        fprintf(fid, '\r\n');
        printflag = 0;
    end
    for j = 1:length(unique_stim_type)
        buff = sprintf('%s\t %5.1f\t %5.2f\t %5.3f\t %5.2f\t %5.2f\t %5.2f\t %5.1f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.3f\t %5.3f\t %10.8f\t %5.3f\t %6.3f\t %10.3f\t %10.3f\t %5.2f\t %5.2f\t', ...
            FILE, data.neuron_params(PREFERRED_ORIENTATION, 1), data.neuron_params(PREFERRED_SPATIAL_FREQ, 1), data.neuron_params(PREFERRED_TEMPORAL_FREQ, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
            null_rate, DMI, p_value, base_rate, amplitude, pref_tempfreq, max_rate, width, stats1(1), stats1(3), stats2(1), stats2(3), DSI, chi2, chiP);
            %unique_stim_type(j), pmax(j).x, pmax(j).y, pmin(j).x, pmin(j).y, avg_resp(j), null_rate, DMI(j), DTI(j), p_value(j), DDI(j), var_term(j), MSgroup(j), MSerror(j), max(unique_tf), min(unique_tf) );
        fprintf(fid, '%s', buff);
        fprintf(fid, '\r\n');
    end
    fclose(fid);
    %------------------------------------------------------------------------

    return;

%OLD CODE:
%
% % Calculate spontaneous rates before looping through so can calculate DTI
% null_resp = data.spike_rates(SpikeChan, null_trials & select_trials);
% null_rate = mean(null_resp);
% 
% 
% 
% for i=1:length(unique_stim_type)	%for each different stim_type value, plot a separate tuning curve
%     stim_type_select = logical( (stim_type == unique_stim_type(i)) );
%     
%     plot_x = tf(stim_type_select & ~null_trials & select_trials);
%     plot_y = spike_rates(stim_type_select & ~null_trials & select_trials); 
%     
%     %compute the Discrimination Index
%     [DDI(i), var_term(i)] = Compute_DDI(plot_x, plot_y);
%           
%     %NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
%     figure(2);
%     hold on;
%     [px, py, perr, pmax(i), pmin(i)] = PlotTuningCurve(plot_x', plot_y', symbols{i}, lines{i}, 1, 1);
%         
%     %Compute DTI from spline fit
%     DTI(i) = 1 - (pmin(i).y - null_rate)/(pmax(i).y - null_rate);
% 
%     %Calculate modulation index using sqrt raw responses and subtracting spontaneous
%     DMI(i) = Compute_ModIndex(plot_x, plot_y, null_resp)
%         
%     %store mean rates for output
%     mean_rates(i, : ) = py';
%     
%     hold on;
%     errorbar(px, py, perr, perr, symbols{i});
%     hold off;
%     
%     % do the ANOVA on sqrt(firing rate) a la Prince et al.
%     [p_value(i), MSgroup(i), MSerror(i)] = spk_anova_F(sqrt(plot_y), plot_x, px);
%     avg_resp(i) = mean(plot_y);      
%     
% end
% 
% %now, get the firing rate for NULL condition trials and add spontaneous rate to plot
% null_x = [min(px) max(px)];
% null_y = [null_rate null_rate];
% hold on;
% plot(null_x, null_y, 'k--');
% hold off;
% 
% yl = YLim;
% YLim([0 yl(2)]);	% set the lower limit of the Y axis to zero
% XLabel('Grating Temporal Frequency(cyc/sec)');
% YLabel('Response (spikes/sec)');
% 
% %now, print out some useful information in the upper subplot
% subplot(2, 1, 1);
% PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
% 
% %now print out useful values 
% % pmax, pmin, py, 
% PrintGratingData(p_value, avg_resp, pmax, pmin, px, null_rate, unique_stim_type, PATH, FILE, DDI, DTI);
% 
% output = 1;
% if (output == 1)
%     
%     %------------------------------------------------------------------------
%     %write out all relevant parameters to a cumulative text file, GCD 8/08/01
%     %write out one line for each stimu_type for each neuron.
%     outfile = [BASE_PATH 'ProtocolSpecific\GratingTemporalFreq\GratingTF_Summary.dat'];
%     printflag = 0;
%     if (exist(outfile, 'file') == 0)    %file does not yet exist
%         printflag = 1;
%     end
%     fid = fopen(outfile, 'a');
%     if (printflag)
%         fprintf(fid, 'FILE\t\t\t PrfOR\t PrfSF\t PrfTF\t RFXctr\t RFYctr\t RFDiam\t GrType\t MaxLoc\t MaxRsp\t MinLoc\t MinRsp\t AvgRsp\t Spont\t ModIx\t DTI  \t AnovaP\t DDI  \t VarTerm\t MSgrp\t MSerror\t MaxX \t MinX \t');
%         fprintf(fid, '\r\n');
%         printflag = 0;
%     end
%     for j = 1:length(unique_stim_type)
%         buff = sprintf('%s\t %5.1f\t %5.2f\t %5.3f\t %5.2f\t %5.2f\t %5.2f\t %5.1f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.3f\t %5.3f\t %10.8f\t %5.3f\t %6.3f\t %10.3f\t %10.3f\t %5.2f\t %5.2f\t', ...
%             FILE, data.neuron_params(PREFERRED_ORIENTATION, 1), data.neuron_params(PREFERRED_SPATIAL_FREQ, 1), data.neuron_params(PREFERRED_TEMPORAL_FREQ, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
%             unique_stim_type(j), pmax(j).x, pmax(j).y, pmin(j).x, pmin(j).y, avg_resp(j), null_rate, DMI(j), DTI(j), p_value(j), DDI(j), var_term(j), MSgroup(j), MSerror(j), max(unique_tf), min(unique_tf) );
%         fprintf(fid, '%s', buff);
%         fprintf(fid, '\r\n');
%     end
%     fclose(fid);
%     %------------------------------------------------------------------------
%     
% end  %if (output == 1)
% 
% return;