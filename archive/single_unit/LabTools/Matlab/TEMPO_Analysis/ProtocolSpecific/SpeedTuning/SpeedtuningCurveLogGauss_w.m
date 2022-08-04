%-----------------------------------------------------------------------------------------------------------------------
%-- SpeedTuningCurveLogGauss.m -- Plots a speed tuning curve.
%--	GCD, 1/26/00
%-----------------------------------------------------------------------------------------------------------------------
function SpeedTuningCurveLogGauss(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

Path_Defs;
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01

% not implemented yet to select output
output = 1;

%get the column of values of speeds in the dots_params matrix
speed = data.dots_params(DOTS_SPEED,:,PATCH1);
%speed'

%now, get the firing rates for all the trials 
spike_rates = data.spike_rates(SpikeChan, :);

%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (speed == data.one_time_params(NULL_VALUE)) );

%now, remove trials from speed and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(speed);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

plot_x = speed(~null_trials & select_trials);
plot_y = spike_rates(~null_trials & select_trials);

figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'Speed Tuning Curve');
subplot(2, 1, 2);

%NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
[px, py, perr] = PlotTuningCurve(plot_x', plot_y', 'ko', 'k-', 1, 0);
errorbar(px, py, perr, perr, 'ko');

% Compute a speed discrimination index analogous to the DDI
[SDI, var_term] = Compute_DDI(plot_x, plot_y);

% Now fit data to loggauss function
means = [px py];
raw = [plot_x' plot_y'];
[pars] = loggaussfit(means, raw);
x_interp = (px(1)): 0.05 : (px(length(px)));
y_interp = loggaussfunc(x_interp, pars);
fit_err = loggauss_err(pars);

%Compute R^2 of the fit for both means and raw values
y_fit = loggaussfunc(px, pars);
y_fit(y_fit < 0) = 0;
%add a column of ones to yfit to make regress happy
y_fit = [ones(length(y_fit),1) y_fit];
[b, bint, r, rint, stats1] = regress(py, y_fit);

y_fit_raw = loggaussfunc(plot_x', pars);
y_fit_raw(y_fit_raw < 0) = 0;
y_fit_raw = [ones(length(y_fit_raw),1) y_fit_raw];
[b, bint, r, rint, stats2] = regress(plot_y', y_fit_raw);

% Do chi-square goodness of fit test
[chi2, chiP] = Chi2_Test(plot_x, plot_y, 'loggaussfunc', pars, length(pars));

[pmax.y max_i] = max(y_interp);
pmax.x = x_interp(max_i);
[pmin.y min_i] = min(y_interp);
pmin.x = x_interp(min_i);

% Now plot fitted curve
hold on;
plot(x_interp, y_interp, 'k-');
hold off;   

%now, get the firing rate for NULL condition trials and add spontaneous rate to plot
null_x = [min(px) max(px)];
null_resp = data.spike_rates(SpikeChan, null_trials & select_trials);
null_rate = mean(null_resp);
null_y = [null_rate null_rate];
hold on;
plot(null_x, null_y, 'k--');
hold off;

%yl = YLim;
YLim([0 yl(2)]);	% set the lower limit of the Y axis to zero
XLabel('Speed of Motion (deg/sec)');
YLabel('Response (spikes/sec)');

%now, print out some useful information in the upper subplot
%subplot(2, 1, 1);
%PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

% first calculate some metrics and stats
p_value = spk_anova(plot_y, plot_x, px);
avg_resp = mean(plot_y);

%now, print out speed tuning curve specific stuff
PrintSpeedDataLogGauss(plot_x, plot_y, pmax, pmin, px, null_rate, pars, fit_err, p_value, avg_resp, stats1, stats2,SDI,chi2,chiP);

%Calculate modulation index using sqrt raw responses and subtracting spontaneous
SMI = Compute_ModIndex(plot_x, plot_y, null_resp);

%now output parameters about the speed tuning curve
    if output == 1
    i = size(PATH,2) - 1;
    while PATH(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
        i = i - 1;
    end   
    PATHOUT = [PATH(1:i) 'Analysis\Tuning\'];
    %PATHOUT = 'C:\MATLAB6p1\work\speedtuning\';
    i = size(FILE,2) - 1;
    while FILE(i) ~='.'
        i = i - 1;
    end
    FILEOUT = [FILE(1:i) 'logspd'];
    
    fileid = [PATHOUT FILEOUT];
    fwriteid = fopen(fileid, 'w');
    
    %fprintf(fwriteid, ' Max Speed	Max Resp Min Speed	Min Resp	Ave Resp	Spont Resp	Curve ANOVA	Fit q(1)	Fit q(2)	Fit q(3)	Fit q(4)	Fit q(5)	Pref Dir	Pref Speed	Mapped Pref H Disp	RF X-Ctr	RF Y-Ctr	RF Diam\n');  
    
    fprintf(fwriteid, '%.6f	%.6f	%.6f	%.6f	%.6f	%.6f	%.6f	%.6f	%.6f	%.6f	%.6f	%.6f	%.2f	%.2f	%.2f	%.2f	%.2f	%.2f\n', pmax.x, pmax.y, pmin.x, pmin.y, avg_resp, null_rate, p_value, pars(1), pars(2), pars(3), pars(4), pars(5), data.one_time_params(PREFERRED_DIRECTION), data.one_time_params(PREFERRED_SPEED), data.one_time_params(PREFERRED_HDISP), data.one_time_params(RF_XCTR), data.one_time_params(RF_YCTR), data.one_time_params(RF_DIAMETER));
    
    fclose(fwriteid);
    
    %also write out data in form suitable for plotting tuning curve with Origin.
    FILEOUT2 = [FILE(1:i) 'logspd_curv'];
    fileid = [PATHOUT FILEOUT2];
    proffid = fopen(fileid, 'w');
    fprintf(proffid,'SpdIn\tFit\tSpeed\tAvgResp\tStdErr\tSpeed2\tSpon\n');
    for i=1:length(x_interp)
        fprintf(proffid,'%6.2f\t%6.2f\t', x_interp(i), y_interp(i));
        if (i <= length(px))
            fprintf(proffid,'%6.2f\t%6.2f\t%6.3f\t', px(i), py(i), perr(i));
        else
            fprintf(proffid,'\t\t\t');
        end
        if (i <= 2)
            fprintf(proffid,'%6.2f\t%6.2f\n',null_x(i),null_y(i));
        else
            fprintf(proffid,'\t\n');
        end
    end
    fclose(proffid);
    
    %---------------------------------------------------------------------------------------
    %ALso, write out some summary data to a cumulative summary file
    if (px(1) == 0)
        zero_resp = py(1)
    else
        zero_resp = NaN;
    end
    buff = sprintf('%s\t %6.1f\t %6.2f\t %6.3f\t %6.2f\t %6.2f\t %6.2f\t %5.1f\t %6.2f\t %6.3f\t %6.2f\t %6.3f\t %6.2f\t %6.2f\t %6.3f\t %10.8f\t %6.3f\t %10.8f\t %6.3f\t %10.8f\t %6.4f\t %6.3f\t %6.2f\t %6.2f\t %10.7f\t %10.7f\t %10.7f\t %8.5f\t %10.8f\t', ...
        FILE, data.neuron_params(PREFERRED_DIRECTION, 1), data.neuron_params(PREFERRED_SPEED, 1), data.neuron_params(PREFERRED_HDISP, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
        0, zero_resp, pmax.x, pmax.y, pmin.x, pmin.y, null_rate, SMI, p_value, stats1(1), stats1(3), stats2(1), stats2(3), SDI, var_term, pars, chi2, chiP);
    outfile = [BASE_PATH 'ProtocolSpecific\SpeedTuning\LogTuningSummary.dat'];
    printflag = 0;
    if (exist(outfile, 'file') == 0)    %file does not yet exist
        printflag = 1;
    end
    fid = fopen(outfile, 'a');
    if (printflag)
        fprintf(fid, 'FILE\t\t PrDir\t PrSpd\t PrHDsp\t RFX\t RFY\t RFDiam\t ZerSpd\t ZerRsp\t MaxSpd\t MaxRsp\t MinSpd\t MinRsp\t Spont\t SpdMI\t AnovaP\t\t Rmeans\t Pmeans\t\t Rraw\t Praw\t\t SDI\t VarTrm\t R0\t K\t vn\t\t sigma\t\t v0\t\t Chi2\t\t ChiP\t\t');
        fprintf(fid, '\r\n');
    end
    fprintf(fid, '%s', buff);
    fprintf(fid, '\r\n');
    fclose(fid);
    %---------------------------------------------------------------------------------------

    %output a cumulative file of the LogGauss function parameters
    outfile2 = [BASE_PATH 'ProtocolSpecific\SpeedTuning\LogParams.dat'];
    if (exist(outfile2, 'file') == 0)    %file does not yet exist
        printflag = 1;
    end
    fsummid = fopen(outfile2, 'a');
    if (printflag)
        fprintf(fsummid, 'FILE\t\t q(1)\t q(2)\t q(3)\t q(4)\t q(5)\t spont\t');
        fprintf(fsummid, '\r\n');
    end
    fprintf(fsummid, '%s\t %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f', FILE, pars(1), pars(2), pars(3), pars(4), pars(5), null_rate);
    fprintf(fsummid, '\r\n');
    fclose(fsummid);
    
end

return;