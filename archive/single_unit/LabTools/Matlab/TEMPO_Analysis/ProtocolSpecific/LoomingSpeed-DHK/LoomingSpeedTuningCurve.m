%-----------------------------------------------------------------------------------------------------------------------
%-- SpeedTuningCurve.m -- Plots a speed tuning curve.
%--	GCD, 1/26/00
%-----------------------------------------------------------------------------------------------------------------------
function LoomingSpeedTuningCurve(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

Path_Defs;
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01

% not implemented yet to select output
output = 1;

%get the column of values of speeds in the dots_params matrix
speed = data.dots_params(DOTS_AP_VEL,:,PATCH1);

%speed'

%now, get the firing rates for all the trials
spike_rates = data.spike_rates(SpikeChan, :);

%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (speed == data.one_time_params(NULL_VALUE)) );
unique_speed = munique(speed(~null_trials)');

%now, remove trials from speed and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(speed);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );



plot_x = speed(~null_trials & select_trials);
plot_y = spike_rates(~null_trials & select_trials);

figure;
set(gcf,'PaperPosition', [.51 .51 20 27.5], 'Position', [250 50 500 573], 'Name', 'Speed Tuning Curve');
subplot(2, 1, 2);
hold on;

%NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
[px, py, perr] = PlotTuningCurve(plot_x', plot_y', 'ko', 'k-', 1, 0);
errorbar(px, py, perr, perr, 'ko');
unique_speeds = px;

% Compute a speed discrimination index analogous to the DDI
[SDI, var_term] = Compute_DDI(plot_x, plot_y);

%Now, fit the data with a log gaussian curve and plot this as well
means = [px py];
raw = [plot_x' plot_y'];
[pars] = speed_loggaussfit(means, raw) %last arg=0 => only allow positive going fit
x_interp = (px(1)):0.01:(px(length(px)));
y_interp = speed_loggaussfunc(x_interp, pars);
plot(x_interp, y_interp, 'k-');


% Now plot fitted curve
hold on;
plot(x_interp, y_interp, 'k-');


yl=ylim;
ylim([0 yl(2)]); %set lower limit of the Y axis to zero
xlabel('Speed of Motion (deg/sec)');
ylabel('Response (spikes/sec)');

%Compute R^2 for both means and raw values
y_fit = speed_loggaussfunc(px, pars);
y_fit(y_fit<0) = 0;
%add a column of ones to yfit to make regress happy
y_fit = [ones(length(y_fit), 1) y_fit];
[b, bint, r, rint, stats1] = regress(py, y_fit);

y_fit_raw = speed_loggaussfunc(plot_x', pars);
y_fit_raw(y_fit_raw<0) = 0;
y_fit_raw = [ones(length(y_fit_raw), 1) y_fit_raw];
[b, bint, r, rint, stats2] = regress(plot_y', y_fit_raw);

%Calculate some metrics and stats to print in the plot
p_value = spk_anova(plot_y, plot_x, unique_speeds);
base_rate = pars(1);
amplitude = pars(2);
peak = pars(3);
st_dev = pars(4);
log_offset = pars(5);
max_rate = base_rate + amplitude;

%now, get the firing rate for NULL condition trials and add spontaneous rate to plot
null_x = [min(px) max(px)];
null_resp = data.spike_rates(SpikeChan, null_trials & select_trials);
null_rate = mean(null_resp);
null_y = [null_rate null_rate];
hold on;
plot(null_x, null_y, 'k--');
hold off;

SMI = 1 - (base_rate - null_rate)/(max_rate - null_rate);

%Do chi-square goodness of fit test
[chi2, chiP] = Chi2_Test(plot_x, plot_y, 'speed_loggaussfunc', pars, length(pars));

[pmax.y max_i] = max(y_interp);
pmax.x = x_interp(max_i);
[pmin.y min_i] = min(y_interp);
pmin.x = x_interp(min_i);

if p_value < 0.05 %& stats1(1) > 0.7;
    index = max_i;
    ysearch = pmax.y;
    HWHMR = 0;
    while (ysearch > (pmin.y + (pmax.y-pmin.y)/2))
        ysearch = y_interp(index);
        if (index == length(y_interp))
            HWHMR = 9999;
            break;
        else
            index = index + 1;
        end
    end
    if HWHMR ~= 9999;
        HWHMR = x_interp(index) - pmax.x;
    end


    index = max_i;
    ysearch = pmax.y;
    HWHML = 0;
    while (ysearch > (pmin.y + (pmax.y-pmin.y)/2))
        ysearch = y_interp(index);
        if (index == 1)
            HWHML = 9999;
            break;
        else
            index = index - 1;
        end
    end
    if HWHML ~= 9999;
        HWHML = pmax.x - x_interp(index);
    end
else
    HWHML = 9999;
    HWHMR = 9999;
end


%now, print out some useful information in the upper subplot
subplot(2, 1, 1);
PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

% first calculate some metrics and stats
p_value = spk_anova(plot_y, plot_x, px);
avg_resp = mean(plot_y);

%now, print out speed tuning curve specific stuff
PrintSpeedData(p_value, base_rate, amplitude, peak, st_dev, log_offset, max_rate, stats1, stats2, chi2, chiP, SDI, HWHMR, HWHML);

%Calculate modulation index using sqrt raw responses and subtracting spontaneous
SMI = Compute_ModIndex(plot_x, plot_y, null_resp);

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
        
    %------------------------------------------------------------------------
    %write out all relevant parameters to a cumulative text file, GCD 8/08/01
    %write out one line for each stimu_type for each neuron.
    outfile = [BASE_PATH 'ProtocolSpecific\LoomingSpeed-DHK\LoomingSpeed_Summary.dat'];
    printflag = 0;
    if (exist(outfile, 'file') == 0)    %file does not yet exist
        printflag = 1;
    end
    fid = fopen(outfile, 'a');
    if (printflag)
        fprintf(fid, 'FILE\t\t\t PrfOR\t PrfSF\t PrfTF\t RFXctr\t RFYctr\t RFDiam\t     MaxLoc\t MaxRsp\t MinLoc\t MinRsp\t Spont\t SMI\t\t AnovaP\t SDI\t\t VarTerm\t MaxX \t MinX \t HWHML\t HWHMR\t RsqMn\t PMean\t RsqRaw\t Praw\t\t Chi2\t\t ChiP\t\t Brat\t\t Ampl\t\t Peak\t\t SD\t\t Offst\t\t');
        fprintf(fid, '\r\n');
        printflag = 0;
    end

    buff = sprintf('%s\t %5.1f\t %5.2f\t %5.3f\t %5.2f\t %5.2f\t %5.2f\t    %5.4f\t %5.4f\t %5.4f\t %5.4f\t %5.4f\t %5.4f\t %5.4f\t %5.4f\t %5.4f\t %5.4f\t %5.4f\t %5.3f\t %5.3f\t %5.4f\t %5.4f\t %5.4f\t %5.4f\t %5.4f\t %5.4f\t %5.4f\t %5.4f\t %5.4f\t %5.4f\t %5.4f\t', ...
        FILE, data.neuron_params(PREFERRED_ORIENTATION, 1), data.neuron_params(PREFERRED_SPATIAL_FREQ, 1), data.neuron_params(PREFERRED_TEMPORAL_FREQ, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
        pmax.x, pmax.y, pmin.x, pmin.y, null_rate, SMI, p_value, SDI, var_term, max(unique_speed), min(unique_speed), HWHML, HWHMR, stats1(1), stats1(3), stats2(1), stats2(3), chi2, chiP, base_rate, amplitude, peak, st_dev, log_offset );
    fprintf(fid, '%s', buff);
    fprintf(fid, '\r\n');

    fclose(fid);
    %------------------------------------------------------------------------

end  %if (output == 1)

return;