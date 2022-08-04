%-----------------------------------------------------------------------------------------------------------------------
%-- OrientationTuningCurve.m -- Plots an orientation tuning curve,for sine
%and square waves
%--	GCD, 1/12/06
%-----------------------------------------------------------------------------------------------------------------------

%%%% 4/27/16 AM added conditional when writing outfile to write into the correct
%%%% directory if on a specific computer; functionality should be
%%%% unchanged if using on any lab computer
%%%% 7/30/16 AM edited formatting of output file columns 
% 3/14/16 AM fixed capitalization to make compatible with new Matlab
%   versions
%%%% 8/23/16 AM added optional 'skip_trials' argin - these trials will not be analyzed;
%%%%    trials not included within BegTrial-EndTrial will still be skipped

function OrientationTuningCurve(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, skip_trials)

TEMPO_Defs;		%needed for defines like IN_T1_WIN_CD
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01
Path_Defs;

symbols = {'ko' 'k*' 'go' 'mo' 'b*' 'r*' 'g*' 'c*'};
lines = {'k-' 'k--' 'g-' 'm-' 'b--' 'r--' 'g--' 'c--'};


%get the column of values of orientations in gratings_params matrix
orientation = data.gratings_params(GRAT_ORIENTATION,:,PATCH1);
%orientation'

%get indices of any NULL conditions (for measuring spontaneous activity)
null_trials = logical( (orientation == data.one_time_params(NULL_VALUE)) );

%get the column of stimulus type values
stim_type = data.gratings_params(GRAT_TYPE,:,PATCH1);
unique_stim_type = munique(stim_type(~null_trials)');

%now, get the firing rates for all the trials 
spike_rates = data.spike_rates(SpikeChan, :);

unique_orientation = munique(orientation(~null_trials)');

%now, remove trials from that do not fall between BegTrial and EndTrial
trials = 1:length(orientation);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

%%%% AM added 8/29/16
if exist('skip_trials','var')
    select_trials(skip_trials) = false;
end

% Calculate spontaneous rates before looping through so can calculate DTI
null_resp = data.spike_rates(SpikeChan, null_trials & select_trials);
null_rate = mean(null_resp);

figure;
set(gcf,'PaperPosition', [.51 .51 20 27.5], 'Position', [450 50 500 573], 'Name', 'Grating Orientation Tuning Curve');
subplot(2, 1, 2);
hold on

for i=1:length(unique_stim_type)	%for each different stim_type value, plot a separate tuning curve
    stim_type_select = logical( (stim_type == unique_stim_type(i)) );
    
    plot_x = orientation(stim_type_select & ~null_trials & select_trials);
    plot_y = spike_rates(stim_type_select & ~null_trials & select_trials); 
    
    %compute the Discrimination Index
    [DDI(i), var_term(i)] = Compute_DDI(plot_x, plot_y);
          
    %NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
    hold on;
    [px, py, perr] = PlotTuningCurve(plot_x', plot_y', symbols{i}, lines{i}, 1, 0);
    
    %Plot the data
    hold on;
    plot(plot_x, plot_y, 'k.');
    errorbar(px, py, perr, perr, symbols{i});
    hold on;
    
    %Now, fit the data with a sum of two Von Mises functions and plot this as
    %well
    means = [px py];
    raw = [plot_x' plot_y'];
    [pars] = von_mises_fit(means, raw)
    x_interp = [(px(1)):1:(px(length(px)))];
    y_interp = von_mises_func(x_interp, pars);
    plot(x_interp, y_interp, 'k-');
        
    %Compute DTI from spline fit
    %DTI(i) = 1 - (pmin(i).y - null_rate)/(pmax(i).y - null_rate);

    %Calculate modulation index using sqrt raw responses and subtracting spontaneous
    DMI(i) = Compute_ModIndex(plot_x, plot_y, null_resp);
        
    %store mean rates for output
    mean_rates(i, : ) = py';
    
    
    % do the ANOVA on sqrt(firing rate) a la Prince et al.
    [p_value(i), MSgroup(i), MSerror(i)] = spk_anova_F(sqrt(plot_y), plot_x, px);
    avg_resp(i) = mean(plot_y);      
    
    [pmax{i}.y max_i] = max(y_interp);
    pmax{i}.x = x_interp(max_i);
    [pmin{i}.y min_i] = min(y_interp);
    pmin{i}.x = x_interp(min_i);
    
    %Calculate the half-width at half height of the highest peak, and
    %return the location of the highest peak
    [max_y{i}, max_index_y{i}] = max(y_interp);
    [min_y{i}, min_index_y{i}] = min(y_interp);
    max_location(i) = x_interp(max_index_y{i});
    index = max_index_y{i};
    ysearch = max_y{i};
    check = 0;
    if p_value(i)<0.01;
        
        %Only find the HWHM and peak location if the cell shows significant
        %orientation tuning (p<0.01), if the tuning is not significant
        %return 9999
        while (ysearch > (min_y{i} + (max_y{i}-min_y{i})/2))
            ysearch = y_interp(index);
            if (index >= length(y_interp))
                HWHM(i) = 9999;
                disp('No HWHM')
                break
            else
                index = index + 1;
            end
        end
        if (check ~= 1)
            HWHM(i) = x_interp(index) - max_location(i);
        end
        
    else
        max_location(i) = 9999;
        HWHM(i) = 9999;
    end
    
    
end
p_value
max_location
HWHM

%now, get the firing rate for NULL condition trials and add spontaneous rate to plot
null_x = [min(px) max(px)];
null_y = [null_rate null_rate];
hold on;
plot(null_x, null_y, 'k--');


yl = ylim;% 3/14/16 AM changed 'YLim' to 'ylim'
ylim([0 yl(2)]);	% set the lower limit of the Y axis to zero  % 3/14/16 AM changed 'YLim' to 'ylim'
xlabel('Grating Orientation(deg)'); % 3/14/16 AM changed 'XLabel' to 'xlabel'
ylabel('Response (spikes/sec)'); % 3/14/16 AM changed 'YLabel' to 'ylabel'

%now, print out some useful information in the upper subplot
subplot(2, 1, 1);
PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
% 
% %now print out useful values 
% % pmax, pmin, py, 
PrintGratingOrientationData(p_value, avg_resp, pmax, pmin, px, null_rate, unique_stim_type, PATH, FILE, DDI, HWHM);

output = 1;
if (output == 1)
    
    %------------------------------------------------------------------------
    %write out all relevant parameters to a cumulative text file, GCD 8/08/01
    %write out one line for each stimu_type for each neuron.
    
    %%%% 4/27/16 AM added conditional below to write into the correct
    %%%% directory if on a specific computer; if on any other lab computer
    %%%% such as 'BURK-LV', the outfile will be written to the usual directory
    if strcmp(getenv('COMPUTERNAME'),'ANDREW') %%% Andrew's MSI laptop name
        outfile = 'C:\Users\AM\Documents\Lab\recordings\GratingOrientationSummary.dat';
    else
        outfile = [BASE_PATH 'ProtocolSpecific\GratingOrientation\GratingOrientationSummary.dat'];
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
        %%% original line: fprintf(fid, 'FILE\t\t\t PrfOR\t PrfSF\t PrfTF\t RFXctr\t RFYctr\t RFDiam\t GrType\t MaxLoc\t MaxRsp\t MinLoc\t MinRsp\t AvgRsp\t Spont\t ModIx\t AnovaP\t DDI  \t VarTerm\t MSgrp\t MSerror\t MaxX \t MinX \t MaxLoc\t HWHM\t');
        fprintf(fid, 'FILE\t PrfOR\t PrfSF\t PrfTF\t RFXctr\t RFYctr\t RFDiam\t GrType\t MaxLoc\t MaxRsp\t MinLoc\t MinRsp\t AvgRsp\t Spont\t ModIx\t AnovaP\t DDI  \t VarTerm\t MSgrp\t MSerror\t MaxX \t MinX \t MaxLoc\t HWHM\t');
        fprintf(fid, '\r\n');
        printflag = 0;
    end
    
    for j = 1:length(unique_stim_type)
        buff = sprintf('%s\t %5.1f\t %5.2f\t %5.3f\t %5.2f\t %5.2f\t %5.2f\t %5.1f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.3f\t %10.8f\t %5.3f\t %6.3f\t %10.3f\t %10.3f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t', ...
            FILE, data.neuron_params(PREFERRED_ORIENTATION, 1), data.neuron_params(PREFERRED_SPATIAL_FREQ, 1), data.neuron_params(PREFERRED_TEMPORAL_FREQ, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
            unique_stim_type(j), pmax{j}.x, pmax{j}.y, pmin{j}.x, pmin{j}.y, avg_resp(j), null_rate, DMI(j), p_value(j), DDI(j), var_term(j), MSgroup(j), MSerror(j), max(unique_orientation), min(unique_orientation), max_location(j), HWHM(j));
        fprintf(fid, '%s', buff);
        fprintf(fid, '\r\n');
    end
    fclose(fid);
    %------------------------------------------------------------------------
    
end  %if (output == 1)

return;