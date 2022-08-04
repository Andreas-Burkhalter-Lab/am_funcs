%-----------------------------------------------------------------------------------------------------------------------
%-- DirectionTuningCurve.m -- Plots a direction tuning curve and computes/plots a Gaussian fit to this curve
%--	GCD, 1/23/00
%-----------------------------------------------------------------------------------------------------------------------

%%%% 4/27/16 AM added conditional when writing outfile to write into the correct
%%%% directory if on a specific computer; functionality should be
%%%% unchanged if using on any lab computer
%%% 7/29/16 AM edited formatting of output file columns; note: using
%%%% Matlab's 'Import Data' option to view the output
%%%% 'DirectionTuningSummary.dat' file may appear to misalign columns; use 
%%%% tdfread for importing or view in excel for proper alignement
%%%% 8/23/16 AM added optional 'skip_trials' argin - these trials will not be analyzed;
%%%%    trials not included within BegTrial-EndTrial will still be skipped

function DirectionTuningCurve(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, skip_trials)

% not implemented yet to select output
output = 1;

Path_Defs;
ProtocolDefs; %contains protocol specific keywords - 1/4/01 BJP

%get the column of values of directions in the dots_params matrix
direction = data.dots_params(DOTS_DIREC,:,PATCH1);

%now, get the firing rates for all the trials 
spike_rates = data.spike_rates(SpikeChan, :);

%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (direction == data.one_time_params(NULL_VALUE)) );

%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(direction);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

%%%% AM added 8/29/16
if exist('skip_trials','var')
    select_trials(skip_trials) = false;
end

plot_x = direction(~null_trials & select_trials);
plot_y = spike_rates(~null_trials & select_trials);

figure;
set(gcf,'PaperPosition', [.51 .51 20 27.5], 'Position', [150 100 500 573], 'Name', 'Direction Tuning Curve');
subplot(2, 1, 2);

%NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
[px, py, perr] = PlotTuningCurve(plot_x', plot_y', 'ko', 'k-', 1, 0);	%note: last arg=0 means just get output, no plot

%keep a copy of the original data before shifting below
px_orig = px;
py_orig = py;

unique_dirs = px; % store category groups for ANOVAs

%now, shift the px, py, and perr vectors such that the peak of tuning curve is in middle of axis range
% now, we need to shift these vectors so that the peak response is always in middle of vector
ctr_indx = round(length(px)/2 - rem(length(px),2)) + 1;
[max_val max_indx] = max(py);
shift = max_indx - ctr_indx;
if (shift > 0)
    px = [ px(shift+1 : length(px)) ; px(1 : shift)+360];
    px_orig = [ px_orig(shift+1 : length(px_orig)) ; px_orig(1 : shift)];
    py = [ py(shift+1 : length(py)) ; py(1 : shift)];
    perr = [ perr(shift+1 : length(perr)) ; perr(1 : shift)];
end
if (shift < 0)
    px = [ px(length(px)+shift+1 : length(px))-360 ; px(1 : length(px)+shift)];
    px_orig = [ px_orig(length(px_orig)+shift+1 : length(px_orig)) ; px_orig(1 : length(px_orig)+shift)];
    py = [ py(length(py)+shift+1 : length(py)) ;  py(1 : length(py)+shift)];
    perr = [ perr(length(perr)+shift+1 : length(perr)) ;  perr(1 : length(perr)+shift)];
end

%now apply the shift to the raw data (plot_x and plot_y)
for i=1:length(px_orig)
    select = logical(plot_x == px_orig(i));
    plot_x(select) = px(i);
end

% Since direction is circular, duplicate  the lowest value (px) as px+ 360, and change spike arrays accordingly
px = [px; px(1)+360];
py = [py; py(1)];
perr = [perr; perr(1)];

hold on;
%plot the data, after shifting as necessary above
plot(plot_x, plot_y, 'k.');
errorbar(px, py, perr, perr, 'ko');
hold on;

% Compute a direction discrimination index analogous to the DDI
[DirDI, var_term] = Compute_DDI(plot_x, plot_y);

%now, fit the data with a Gaussian curve and plot this as well
means = [px py];
raw = [plot_x' plot_y'];
[pars] = gaussfit(means, raw, 0);   %last arg: allow positive going fit only
x_interp = (px(1)): 0.5 : (px(length(px)));
y_interp = gaussfunc(x_interp, pars);
plot(x_interp, y_interp, 'k-');

%now, get the firing rate for NULL condition trials and add spontaneous rate to plot
null_x = [min(px) max(px)];
null_resp = data.spike_rates(SpikeChan, null_trials & select_trials);
null_rate = mean(null_resp);
null_y = [null_rate null_rate];
hold on;
plot(null_x, null_y, 'k--');
hold off;

yl = ylim;
ylim([0 yl(2)]);	% set the lower limit of the Y axis to zero
xlabel('Direction of Motion (deg)');
ylabel('Response (spikes/sec)');

%Compute R^2 of the fit for both means and raw values
y_fit = gaussfunc(px, pars);
y_fit(y_fit < 0) = 0;
%add a column of ones to yfit to make regress happy
y_fit = [ones(length(y_fit),1) y_fit];
[b, bint, r, rint, stats1] = regress(py, y_fit);

y_fit_raw = gaussfunc(plot_x', pars);
y_fit_raw(y_fit_raw < 0) = 0;
y_fit_raw = [ones(length(y_fit_raw),1) y_fit_raw];
[b, bint, r, rint, stats2] = regress(plot_y', y_fit_raw);

% Do chi-square goodness of fit test
[chi2, chiP] = Chi2_Test(plot_x, plot_y, 'gaussfunc', pars, length(pars));

%now, print out some useful information in the upper subplot
subplot(2, 1, 1);
PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

% calculate some metrics and stats then print them in plot
pref_dir = pars(3);
p_value = spk_anova(plot_y, plot_x, unique_dirs);
base_rate = pars(1);
amplitude = pars(2);
max_rate = base_rate + amplitude;
width = sqrt(-(log(.5)))*pars(4)*2*sqrt(2);
DSI = 1 - (base_rate - null_rate)/(max_rate - null_rate); 

%Calculate modulation index using sqrt raw responses and subtracting spontaneous
DMI = Compute_ModIndex(plot_x, plot_y, null_resp);

PrintDirectionData(p_value, base_rate, null_rate, amplitude, pref_dir, max_rate, width, DSI, stats1, stats2, DirDI, chi2, chiP); 


%output tuning curve metrics
if (output == 1)
%     i = size(PATH,2) - 1;
%     while PATH(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
%         i = i - 1;
%     end   
%     PATHOUT = [PATH(1:i) 'Analysis\Tuning\'];
%     i = size(FILE,2) - 1;
%     while FILE(i) ~='.'
%         i = i - 1;
%     end
%     FILEOUT = [FILE(1:i) 'dir'];
%     
%     fileid = [PATHOUT FILEOUT];
%     fwriteid = eval(['fopen(fileid, ''w'')']);
%     %fprintf(fwriteid, '%%Base Rate (q1)	Amplitude (q2)	Pref dir (q3)	q4	Width(FWHM)	Max Resp	Spont Resp	DSI	Curve ANOVA	Mapped Pref Dir	Pref Speed	Pref H Disp	RF X-Ctr	RF Y-Ctr	Diam\n');
%     fprintf(fwriteid, '%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f\n', pars(1), pars(2), pars(3), pars(4), width, max_rate, null_rate, DSI, p_value, data.one_time_params(PREFERRED_DIRECTION), data.one_time_params(PREFERRED_SPEED), data.one_time_params(PREFERRED_HDISP), data.one_time_params(RF_XCTR), data.one_time_params(RF_YCTR), data.one_time_params(RF_DIAMETER));
% 
% 	fclose(fwriteid);
% 
%    %---------------------------------------------------------------------------------------
%    %also write out data in form suitable for plotting tuning curve with Origin.
%     FILEOUT2 = [FILE(1:i) 'direc_curv_fit'];
%     fileid = [PATHOUT FILEOUT2];
%     proffid = fopen(fileid, 'w');
%     fprintf(proffid,'DirIn\tFit\tDirec\tAvgResp\tStdErr\tDir2\tSpon\n');
%     for kk=1:length(x_interp)
%         fprintf(proffid,'%6.2f\t%6.2f\t', x_interp(kk), y_interp(kk));
%         if (kk <= length(px))
%             fprintf(proffid,'%6.2f\t%6.2f\t%6.3f\t', px(kk), py(kk), perr(kk));
%         else
%             fprintf(proffid,'\t\t\t');
%         end
%         if (kk <= 2)
%             fprintf(proffid,'%6.2f\t%6.2f\n',null_x(kk),null_y(kk));
%         else
%             fprintf(proffid,'\t\n');
%         end
%     end
%     fclose(proffid);
    
    %---------------------------------------------------------------------------------------
    %ALso, write out summary data to a cumulative summary file
    [pdir360, pdir180, pdir90] = AngleWrap(pref_dir);
    buff = sprintf('%s\t %6.1f\t %6.2f\t %6.3f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.3f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.3f\t %10.8f\t %6.3f\t %10.8f\t %6.3f\t %10.8f\t %6.4f\t %6.3f\t %8.5f\t %10.8f\t', ...
        FILE, data.neuron_params(PREFERRED_DIRECTION, 1), data.neuron_params(PREFERRED_SPEED, 1), data.neuron_params(PREFERRED_HDISP, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
        null_rate, DMI, base_rate, amplitude, pdir360, pdir180, pdir90, width, DSI, p_value, stats1(1), stats1(3), stats2(1), stats2(3), DirDI, var_term, chi2, chiP);
    
    
    %%%% 4/27/16 AM added conditional below to write into the correct
    %%%% directory if on a specific computer; if on any other lab computer
    %%%% such as 'BURK-LV', the outfile will be written to the usual directory
    if strcmp(getenv('COMPUTERNAME'),'ANDREW') %%% Andrew's MSI laptop name
        outfile = 'C:\Users\AM\Documents\Lab\recordings\DirectionTuningSummary.dat';
    else
        outfile = [BASE_PATH 'ProtocolSpecific\DirectionTuning\DirectionTuningSummary.dat'];
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
        %%% original line: fprintf(fid, 'FILE\t\t PrDir\t PrSpd\t PrHDsp\t RFX\t RFY\t RFDiam\t Spont\t DirMI\t BRate\t Ampl\t PD360\t PD180\t PD90\t FWHM\t DSI\t AnovaP\t\t Rmeans\t Pmeans\t\t Rraw\t Praw\t\t DirDI\t VarTrm\t Chi2\t\t ChiP\t\t');
        fprintf(fid, 'FILE\t PrDir\t PrSpd\t PrHDsp\t RFX\t RFY\t RFDiam\t Spont\t DirMI\t BRate\t Ampl\t PD360\t PD180\t PD90\t FWHM\t DSI\t AnovaP\t Rmeans\t Pmeans\t Rraw\t Praw\t DirDI\t VarTrm\t Chi2\t ChiP\t');
        fprintf(fid, '\r\n');
    end
    fprintf(fid, '%s', buff);
    fprintf(fid, '\r\n');
    fclose(fid);
    %---------------------------------------------------------------------------------------
    
    %output a cumulative file of the Gaussian fit parameters
    outfile = [BASE_PATH 'ProtocolSpecific\DirectionTuning\DirectionParams.dat'];
    printflag = 0;
    if (exist(outfile, 'file') == 0)    %file does not yet exist
        printflag = 1;
    end
    fsummid = fopen(outfile, 'a');
    if (printflag)
        fprintf(fsummid, 'FILE\t\t q(1)\t q(2)\t q(3)\t q(4)\t spont\t');
        fprintf(fsummid, '\r\n');
    end
    fprintf(fsummid, '%s\t %7.5f %7.5f %7.5f %7.5f %7.5f', FILE, pars(1), pars(2), pars(3), pars(4), null_rate);
    fprintf(fsummid, '\r\n');
    fclose(fsummid);
    
end

%print(2); % Uncomment for autoprinting.  JWN 081605
%close(2); % Uncomment for autoprinting.

return;