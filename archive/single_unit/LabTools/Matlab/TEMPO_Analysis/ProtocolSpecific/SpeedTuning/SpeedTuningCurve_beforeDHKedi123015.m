%-----------------------------------------------------------------------------------------------------------------------
%-- SpeedTuningCurve.m -- Plots a speed tuning curve.
%--	GCD, 1/26/00
%-----------------------------------------------------------------------------------------------------------------------
function SpeedTuningCurve(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

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
    
    %also write out data in form suitable for plotting tuning curve with Origin.
    FILEOUT2 = [FILE(1:i) 'spd_curv'];
    fileid = [PATHOUT FILEOUT2];
    proffid = fopen(fileid, 'w');
    fprintf(proffid,'SpdIn\tFit\tSpeed\tAvgResp\tStdErr\tSpeed2\tSpon\n');
    for i=1:length(x_interp)
        fprintf(proffid,'%6.3f\t%6.3f\t', x_interp(i), y_interp(i));
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
    
    %------------------------------------------------------------------------
    %write out all relevant parameters to a cumulative text file, GCD 8/08/01
    %write out one line for each stimu_type for each neuron.
    outfile = [BASE_PATH 'ProtocolSpecific\SpeedTuning\Speed_Summary.dat'];
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Old code from when fit was with a gamma function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Now fit data to gamma function
% means = [px py];
% raw = [plot_x' plot_y'];
% [pars] = gammafit(means, raw);
% incr = 0.05;
% %x_interp = (px(1)): incr : (px(length(px))); %KLUGE for now
% x_interp = 0: incr : max(speed);
% y_interp = gammafunc(x_interp, pars);
% fit_err = gamma_err(pars);

% %Compute R^2 of the fit for both means and raw values
% y_fit = gammafunc(px, pars);
% y_fit(y_fit < 0) = 0;
% %add a column of ones to yfit to make regress happy
% y_fit = [ones(length(y_fit),1) y_fit];
% [b, bint, r, rint, stats1] = regress(py, y_fit);

% y_fit_raw = gammafunc(plot_x', pars);
% y_fit_raw(y_fit_raw < 0) = 0;
% y_fit_raw = [ones(length(y_fit_raw),1) y_fit_raw];
% [b, bint, r, rint, stats2] = regress(plot_y', y_fit_raw);

% %Compute the derivative of the best-fit tuning curve
% deriv = diff(y_interp)/incr;

% % Do chi-square goodness of fit test
% [chi2, chiP] = Chi2_Test(plot_x, plot_y, 'gammafunc', pars, length(pars))

%plot the derivative function
% hold on;
% plot(x_interp(2:(length(x_interp)-1)), abs(deriv(2:length(deriv))),
% 'g--');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Old file output code:
% %------------------------------------------------------------
% %now output parameters about the speed tuning curve
% if output == 1
%     i = size(PATH,2) - 1;
%     while PATH(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
%         i = i - 1;
%     end   
%     PATHOUT = [PATH(1:i) 'Analysis\Tuning\'];
%     i = size(FILE,2) - 1;
%     while FILE(i) ~='.'
%         i = i - 1;
%     end
%     FILEOUT = [FILE(1:i) 'spd'];
%     
%     fileid = [PATHOUT FILEOUT]
%     fwriteid = fopen(fileid, 'w')
%     
%     %fprintf(fwriteid, ' Max Speed	Max Resp Min Speed	Min Resp	Ave Resp	Spont Resp	Curve ANOVA	Fit q(1)	Fit q(2)	Fit q(3)	Fit q(4)	Fit q(5)	Pref Dir	Pref Speed	Mapped Pref H Disp	RF X-Ctr	RF Y-Ctr	RF Diam\n');  
%     
%     fprintf(fwriteid, '%.6f	%.6f	%.6f	%.6f	%.6f	%.6f	%.6f	%.6f	%.6f	%.6f	%.6f	%.6f	%.2f	%.2f	%.2f	%.2f	%.2f	%.2f\n', pmax.x, pmax.y, pmin.x, pmin.y, avg_resp, null_rate, p_value, pars(1), pars(2), pars(3), pars(4), pars(5), data.one_time_params(PREFERRED_DIRECTION), data.one_time_params(PREFERRED_SPEED), data.one_time_params(PREFERRED_HDISP), data.one_time_params(RF_XCTR), data.one_time_params(RF_YCTR), data.one_time_params(RF_DIAMETER));
%     
%     fclose(fwriteid);
%     
%     %also write out data in form suitable for plotting tuning curve with Origin.
%     FILEOUT2 = [FILE(1:i) 'spd_curv'];
%     fileid = [PATHOUT FILEOUT2];
%     proffid = fopen(fileid, 'w');
%     fprintf(proffid,'SpdIn\tFit\tSpeed\tAvgResp\tStdErr\tSpeed2\tSpon\n');
%     for i=1:length(x_interp)
%         fprintf(proffid,'%6.2f\t%6.2f\t', x_interp(i), y_interp(i));
%         if (i <= length(px))
%             fprintf(proffid,'%6.2f\t%6.2f\t%6.3f\t', px(i), py(i), perr(i));
%         else
%             fprintf(proffid,'\t\t\t');
%         end
%         if (i <= 2)
%             fprintf(proffid,'%6.2f\t%6.2f\n',null_x(i),null_y(i));
%         else
%             fprintf(proffid,'\t\n');
%         end
%     end
%     fclose(proffid);
%     
%     %---------------------------------------------------------------------------------------
%     %ALso, write out some summary data to a cumulative summary file
%     if (px(1) == 0)
%         zero_resp = py(1)
%     else
%         zero_resp = NaN;
%     end
%     buff = sprintf('%s\t %6.1f\t %6.2f\t %6.3f\t %6.2f\t %6.2f\t %6.2f\t %5.1f\t %6.2f\t %6.3f\t %6.2f\t %6.3f\t %6.2f\t %6.2f\t %6.3f\t %10.8f\t %6.3f\t %10.8f\t %6.3f\t %10.8f\t %6.4f\t %6.3f\t %6.2f\t %6.2f\t %10.7f\t %10.7f\t %10.7f\t %8.5f\t %10.8f\t %6d\t %6d\t %6.4f\t', ...
%         FILE, data.neuron_params(PREFERRED_DIRECTION, 1), data.neuron_params(PREFERRED_SPEED, 1), data.neuron_params(PREFERRED_HDISP, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
%         0, zero_resp, pmax.x, pmax.y, pmin.x, pmin.y, null_rate, SMI, p_value, stats1(1), stats1(3), stats2(1), stats2(3), SDI, var_term, pars, chi2, chiP, length(plot_y), length(pars), fit_err);
%     outfile = [BASE_PATH 'ProtocolSpecific\SpeedTuning\SpeedTuningSummary_Gamma.dat'];
%     printflag = 0;
%     if (exist(outfile, 'file') == 0)    %file does not yet exist
%         printflag = 1;
%     end
%     fid = fopen(outfile, 'a');
%     if (printflag)
%         fprintf(fid, 'FILE\t\t PrDir\t PrSpd\t PrHDsp\t RFX\t RFY\t RFDiam\t ZerSpd\t ZerRsp\t MaxSpd\t MaxRsp\t MinSpd\t MinRsp\t Spont\t SpdMI\t AnovaP\t\t Rmeans\t Pmeans\t\t Rraw\t Praw\t\t SDI\t VarTrm\t R0\t K\t alpha\t\t tau\t\t expon\t\t Chi2\t\t ChiP\t\t Npts\t Npars\t FitErr\t');
%         fprintf(fid, '\r\n');
%     end
%     fprintf(fid, '%s', buff);
%     fprintf(fid, '\r\n');
%     fclose(fid);
%     %---------------------------------------------------------------------------------------
% 
%    %---------------------------------------------------------------------------------------
%     %write out some speed tuning slope data, obtained from raw rates
%     raw_slope = [];
%     for k=1:(length(px)-1)
%         if (px(k+1)< 35.0)
%             raw_slope(k) = abs((py(k+1)-py(k))/(px(k+1)-px(k)));
%         end
%     end
%     
%     buff = sprintf('%8.5f\t', raw_slope);
%     outfile = [BASE_PATH 'ProtocolSpecific\SpeedTuning\SpeedTuningRawSlopes.dat'];
%     fid = fopen(outfile, 'a');
%     fprintf(fid, '%s', buff);
%     fprintf(fid, '\r\n');
%     fclose(fid);
%     %---------------------------------------------------------------------------------------
%     
%     
%     %output a cumulative file of the Gamma function parameters
%     outfile2 = [BASE_PATH 'ProtocolSpecific\SpeedTuning\SpeedParams.dat'];
%     if (exist(outfile2, 'file') == 0)    %file does not yet exist
%         printflag = 1;
%     end
%     fsummid = fopen(outfile2, 'a');
%     if (printflag)
%         fprintf(fsummid, 'FILE\t\t q(1)\t q(2)\t q(3)\t q(4)\t q(5)\t spont\t');
%         fprintf(fsummid, '\r\n');
%     end
%     fprintf(fsummid, '%s\t %7.5f %7.5f %7.5f %7.5f %7.5f %7.5f', FILE, pars(1), pars(2), pars(3), pars(4), pars(5), null_rate);
%     fprintf(fsummid, '\r\n');
%     fclose(fsummid);
%     
% end

