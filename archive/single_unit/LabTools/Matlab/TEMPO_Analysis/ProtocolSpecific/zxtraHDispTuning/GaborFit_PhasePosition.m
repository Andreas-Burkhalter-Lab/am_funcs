%-----------------------------------------------------------------------------------------------------------------------
%-- GaborFit_PhasePosition.m -- Fits a disparity tuning curve with a full Gabor, or with a restricted Gabor where 
%-- either the position or phase is constrained to zero. NOTE: responses are square-rooted in the error functions, to help homogenize
%-- variance (a la Cumming).  A constrained version of Levenberg-Marquardt optimization is used, as
%-- implemented using 'fmincon'.
%--	GCD, 1/26/00
%-----------------------------------------------------------------------------------------------------------------------
function GaborFit_PhasePosition(fit_type, data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01
Path_Defs;

% not implemented yet to select output
output = 0;

symbols = {'ko' 'k*' 'go' 'mo' 'b*' 'r*' 'g*' 'c*'};
lines = {'k-' 'k--' 'g-' 'm-' 'b--' 'r--' 'g--' 'c--'};

%get the column of values of horiz. disparity in the dots_params matrix
hor_disp = data.dots_params(DOTS_HDISP,:,PATCH1);

%get indices of any NULL conditions (for measuring spontaneous activity)
null_trials = logical( (hor_disp == data.one_time_params(NULL_VALUE)) );

%get the column of speed values
speed = data.dots_params(DOTS_SPEED,:,PATCH1);
unique_speed = munique(speed(~null_trials)');

%now, get the firing rates for all the trials 
spike_rates = data.spike_rates(SpikeChan, :);

%get indices of monoc. and uncorrelated controls
control_trials = logical( (hor_disp == LEYE_CONTROL) | (hor_disp == REYE_CONTROL) | (hor_disp == UNCORR_CONTROL) );

%now, remove trials from hor_disp and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(hor_disp);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

% Calculate spontaneous rates before looping through so can calculate DTI
null_rate = mean(data.spike_rates(SpikeChan, null_trials & select_trials));

figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'HDisp Tuning Curve w/ Best Fit');

for i=length(unique_speed):length(unique_speed)	%for each different speed value, plot a separate disparity tuning curve
    %NOTE; right now this is kludged to analyze only the largest (e.g. non-zero) speed.  It may do crazy
    %things otherwise in terms of plotting.  GCD, 1/31/01
    speed_select = logical( (speed == unique_speed(i)) );
    
    plot_x = hor_disp(speed_select & ~null_trials & ~control_trials & select_trials);
    plot_y = spike_rates(speed_select & ~null_trials & ~control_trials & select_trials); 
    
    %NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
    [px, py, perr, pmax(i), pmin(i)] = PlotTuningCurve(plot_x', plot_y', symbols{i}, lines{i}, 1, 0);
    
    p_value(i) = spk_anova(plot_y, plot_x, px);
    
    %now, compute the monoc. and uncorrelated control values
    Leye_trials = logical( (hor_disp == LEYE_CONTROL) );
    Leye_resp = spike_rates(speed_select & select_trials & Leye_trials);
    Reye_trials = logical( (hor_disp == REYE_CONTROL) );
    Reye_resp = spike_rates(speed_select & select_trials & Reye_trials);
    Uncorr_trials = logical( (hor_disp == UNCORR_CONTROL) );
    Uncorr_resp = spike_rates(speed_select & select_trials & Uncorr_trials);
    
    means = [px py];
    raw = [plot_x' plot_y'];
    
    %-------------------------------
    %fit the FULL GABOR FUNCTION
    %-------------------------------
    subplot(4,2,2);
    
    %plot the raw data
    PlotRawData(plot_x', plot_y', symbols{i}, max(px), Leye_resp, Leye_trials, Reye_resp, Reye_trials, Uncorr_resp, Uncorr_trials, null_rate);
    
    fixed_param_flags = zeros(6,1); %by default, all 6 parameters will vary
    fixed_param_values = zeros(6,1); %override these values and flags to fix a parameter
    [pars{i},freq(i)] = gaborfit(means,raw,fixed_param_flags,fixed_param_values);
    
    x_interp = (px(1)): .01 : (px(length(px)));
    y_gabor = gaborfunc(x_interp, pars{i});
    y_gabor(y_gabor < 0) = 0;
    y_gauss =  pars{i}(1) + pars{i}(2)*exp(-0.5*((x_interp - pars{i}(3))/ pars{i}(4)).^2);
    y_sine =  pars{i}(1) + 0.5*pars{i}(2)*cos(2*pi*pars{i}(5)*(x_interp - pars{i}(3))+pars{i}(6) );
    
    % Now plot fitted curves; square the fitted values to counteract sqrt() above
    hold on;
    plot(x_interp, y_gabor, 'k-');
    plot(x_interp, y_gauss, 'k--');
    plot(x_interp, y_sine, 'k:');
    Gabor_SSE(i) = gaborerr(pars{i});
    buff = sprintf('Full Gabor fit (sse = %8.5f)', gaborerr(pars{i}) );
    title(buff);
    hold off;   
    
    y_fit = gaborfunc(px, pars{i});
    y_fit(y_fit < 0) = 0;
    %add a column of ones to yfit to make regress happy
    y_fit = [ones(length(y_fit),1) y_fit];
    [b, bint, r, rint, stats1{i}] = regress(py, y_fit);
    
    y_fit_raw = gaborfunc(plot_x', pars{i});
    y_fit_raw(y_fit_raw < 0) = 0;
    y_fit_raw = [ones(length(y_fit_raw),1) y_fit_raw];
    [b, bint, r, rint, stats2{i}] = regress(plot_y', y_fit_raw);
    
    % do a chi-square goodness of fit test
    [chi2, chiP] = chi2_test(plot_x, plot_y, 'gaborfunc', pars{i}, (length(pars{i})-1) );
    
    %print fit parameters to the screen
    subplot(4,2,1);
    axis([0 100 0 100]);    axis('off');
    xpos = -10; ypos = 100;
    font_size = 7;  bump_size = 8;
    temp = strcat(PATH, FILE);
    temp(temp == '\') = '/';
    % this prevents a stupid error from appearing on the screen
    line = sprintf('File: %s', temp);
    text(xpos-15,ypos, line,'FontSize',font_size+1);		ypos = ypos - 2*bump_size;
    line = sprintf('Base Rate: %8.4f', pars{i}(1));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Amplitude: %8.4f', pars{i}(2));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Gauss Ctr: %8.4f', pars{i}(3));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Gauss  SD: %8.4f', pars{i}(4));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Sine Freq: %8.4f', pars{i}(5));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Sin Phase: %8.4f', pars{i}(6));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('FT Freq: %8.4f', freq(i));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Means: Rsq=%6.3f, P=%8.6f', stats1{i}(1), stats1{i}(3));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Raw:   Rsq=%6.3f, P=%8.6f', stats2{i}(1), stats2{i}(3));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Fit test: chi2=%6.3f, chiP=%8.6f', chi2, chiP);
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    
    % put data into strings for output to file
    names_genl =  sprintf('FILE\t\t PrDir\t PrSpd\t PrHDsp\t RFx\t RFy\t RFdiam\t');
    outstr_genl = sprintf('%s\t %6.1f\t %6.2f\t %6.3f\t %6.2f\t %6.2f\t %6.2f\t', ...
        FILE, data.neuron_params(PREFERRED_DIRECTION, 1), data.neuron_params(PREFERRED_SPEED, 1), data.neuron_params(PREFERRED_HDISP, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1) );
        
    [phase360, phase180, phase90] = AngleWrap(pars{i}(6)*180/pi);
    if (phase360 > 180)
        phase360 = phase360 - 360;  %fold into -180->180 range
    end
    
    names_FullGabor =  sprintf('Rbase1\t Amp1\t Gctr1\t Gwid1\t Sfreq1\t Sphs1\t Sph360\t Sph180\t Sph90\t FTfrq1\t R2mn1\t Pmean1\t\t R2raw1\t Praw1\t\t chisq1\t chiP1\t\t');
    outstr_FullGabor = sprintf('%6.2f\t %6.2f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.1f\t %6.1f\t %6.1f\t %6.3f\t %6.4f\t %10.8f\t %6.4f\t %10.8f\t %6.2f\t %10.8f\t', ...
        pars{i}, phase360, phase180, phase90, freq(i), stats1{i}(1), stats1{i}(3), stats2{i}(1), stats2{i}(3), chi2, chiP);
    
    
    %--------------------------------------
    %fit the POSITION-FIXED GABOR FUNCTION
    %--------------------------------------
    subplot(4,2,4);
    
    %plot the raw data
    PlotRawData(plot_x', plot_y', symbols{i}, max(px), Leye_resp, Leye_trials, Reye_resp, Reye_trials, Uncorr_resp, Uncorr_trials, null_rate);
    
    fixed_param_flags = zeros(6,1); %by default, all 6 parameters will vary
    fixed_param_values = zeros(6,1); %override these values and flags to fix a parameter
    %set the position to be fixed to zero
    fixed_param_flags(3) = 1;
    fixed_param_values(3) = 0;
    [pars{i},freq(i)] = gaborfit(means,raw,fixed_param_flags,fixed_param_values);
    
    x_interp = (px(1)): .01 : (px(length(px)));
    y_gabor = gaborfunc(x_interp, pars{i});
    y_gabor(y_gabor < 0) = 0;
    y_gauss =  pars{i}(1) + pars{i}(2)*exp(-0.5*((x_interp - pars{i}(3))/ pars{i}(4)).^2);
    y_sine =  pars{i}(1) + 0.5*pars{i}(2)*cos(2*pi*pars{i}(5)*(x_interp - pars{i}(3))+pars{i}(6) );
    
    % Now plot fitted curves; square the fitted values to counteract sqrt() above
    hold on;
    plot(x_interp, y_gabor, 'k-');
    plot(x_interp, y_gauss, 'k--');
    plot(x_interp, y_sine, 'k:');
    Gabor_SSE_pos_fixed(i) = gaborerr(pars{i});
    buff = sprintf('Zero Position (sse = %8.5f)', gaborerr(pars{i}) );
    title(buff);
    hold off;   
    
    y_fit = gaborfunc(px, pars{i});
    y_fit(y_fit < 0) = 0;
    %add a column of ones to yfit to make regress happy
    y_fit = [ones(length(y_fit),1) y_fit];
    [b, bint, r, rint, stats1{i}] = regress(py, y_fit);
    
    y_fit_raw = gaborfunc(plot_x', pars{i});
    y_fit_raw(y_fit_raw < 0) = 0;
    y_fit_raw = [ones(length(y_fit_raw),1) y_fit_raw];
    [b, bint, r, rint, stats2{i}] = regress(plot_y', y_fit_raw);
    
    % do a chi-square goodness of fit test
    [chi2_PosFixed, chiP_PosFixed] = chi2_test(plot_x, plot_y, 'gaborfunc', pars{i}, (length(pars{i})-1) );

    %print fit parameters to the screen
    subplot(4,2,3);
    axis([0 100 0 100]);    axis('off');
    xpos = -10; ypos = 100;
    font_size = 7;  bump_size = 8;
    line = sprintf('Base Rate: %8.4f', pars{i}(1));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Amplitude: %8.4f', pars{i}(2));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Gauss Ctr: %8.4f', pars{i}(3));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Gauss  SD: %8.4f', pars{i}(4));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Sine Freq: %8.4f', pars{i}(5));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Sin Phase: %8.4f', pars{i}(6));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('FT Freq: %8.4f', freq(i));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Means: Rsq=%6.3f, P=%8.6f', stats1{i}(1), stats1{i}(3));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Raw:   Rsq=%6.3f, P=%8.6f', stats2{i}(1), stats2{i}(3));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Fit test: chi2=%6.3f, chiP=%8.6f', chi2_PosFixed, chiP_PosFixed);
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;

    [phase360, phase180, phase90] = AngleWrap(pars{i}(6)*180/pi);
    if (phase360 > 180)
        phase360 = phase360 - 360;  %fold into -180->180 range
    end

    names_PosFixedGabor =  sprintf('Rbase2\t Amp2\t Gctr2\t Gwid2\t Sfreq2\t Sphas2\t Sph360\t Sph180\t Sph90\t FTfrq2\t R2mn2\t Pmean2\t\t R2raw2\t Praw2\t\t chisq2\t chiP2\t\t');
    outstr_PosFixedGabor = sprintf('%6.2f\t %6.2f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.1f\t %6.1f\t %6.1f\t %6.3f\t %6.4f\t %10.8f\t %6.4f\t %10.8f\t %6.2f\t %10.8f\t', ...
        pars{i}, phase360, phase180, phase90, freq(i), stats1{i}(1), stats1{i}(3), stats2{i}(1), stats2{i}(3), chi2_PosFixed, chiP_PosFixed);
        
    %--------------------------------------
    %fit the PHASE-FIXED GABOR FUNCTION
    %--------------------------------------
    subplot(4,2,6);
    
    %plot the raw data
    PlotRawData(plot_x', plot_y', symbols{i}, max(px), Leye_resp, Leye_trials, Reye_resp, Reye_trials, Uncorr_resp, Uncorr_trials, null_rate);
    
    fixed_param_flags = zeros(6,1); %by default, all 6 parameters will vary
    fixed_param_values = zeros(6,1); %override these values and flags to fix a parameter
    %set the position to be fixed to zero
    fixed_param_flags(6) = 1;
    fixed_param_values(6) = 0;
    [pars{i},freq(i)] = gaborfit(means,raw,fixed_param_flags,fixed_param_values);
    
    x_interp = (px(1)): .01 : (px(length(px)));
    y_gabor = gaborfunc(x_interp, pars{i});
    y_gabor(y_gabor < 0) = 0;
    y_gauss =  pars{i}(1) + pars{i}(2)*exp(-0.5*((x_interp - pars{i}(3))/ pars{i}(4)).^2);
    y_sine =  pars{i}(1) + 0.5*pars{i}(2)*cos(2*pi*pars{i}(5)*(x_interp - pars{i}(3))+pars{i}(6) );
    
    % Now plot fitted curves; square the fitted values to counteract sqrt() above
    hold on;
    plot(x_interp, y_gabor, 'k-');
    plot(x_interp, y_gauss, 'k--');
    plot(x_interp, y_sine, 'k:');
    Gabor_SSE_phase_fixed(i) = gaborerr(pars{i});
    buff = sprintf('Zero Phase (sse = %8.5f)', gaborerr(pars{i}) );
    title(buff);
    hold off;   
    
    y_fit = gaborfunc(px, pars{i});
    y_fit(y_fit < 0) = 0;
    %add a column of ones to yfit to make regress happy
    y_fit = [ones(length(y_fit),1) y_fit];
    [b, bint, r, rint, stats1{i}] = regress(py, y_fit);
    
    y_fit_raw = gaborfunc(plot_x', pars{i});
    y_fit_raw(y_fit_raw < 0) = 0;
    y_fit_raw = [ones(length(y_fit_raw),1) y_fit_raw];
    [b, bint, r, rint, stats2{i}] = regress(plot_y', y_fit_raw);

    % do a chi-square goodness of fit test
    [chi2_PhaseFixed, chiP_PhaseFixed] = chi2_test(plot_x, plot_y, 'gaborfunc', pars{i}, (length(pars{i})-1) );
    
    %print fit parameters to the screen
    subplot(4,2,5);
    axis([0 100 0 100]);    axis('off');
    xpos = -10; ypos = 100;
    font_size = 7;  bump_size = 8;
    line = sprintf('Base Rate: %8.4f', pars{i}(1));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Amplitude: %8.4f', pars{i}(2));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Gauss Ctr: %8.4f', pars{i}(3));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Gauss  SD: %8.4f', pars{i}(4));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Sine Freq: %8.4f', pars{i}(5));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Sin Phase: %8.4f', pars{i}(6));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('FT Freq: %8.4f', freq(i));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Means: Rsq=%6.3f, P=%8.6f', stats1{i}(1), stats1{i}(3));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Raw:   Rsq=%6.3f, P=%8.6f', stats2{i}(1), stats2{i}(3));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Fit test: chi2=%6.3f, chiP=%8.6f', chi2_PhaseFixed, chiP_PhaseFixed);
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;

    names_PhaseFixedGabor =  sprintf('Rbase3\t Amp3\t Gctr3\t Gwid3\t Sfreq3\t Sphas3\t FTfrq3\t R2mn3\t Pmean3\t\t R2raw3\t Praw3\t\t chisq3\t chiP3\t\t');
    outstr_PhaseFixedGabor = sprintf('%6.2f\t %6.2f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.4f\t %10.8f\t %6.4f\t %10.8f\t %6.2f\t %10.8f\t', ...
        pars{i}, freq(i), stats1{i}(1), stats1{i}(3), stats2{i}(1), stats2{i}(3), chi2_PhaseFixed, chiP_PhaseFixed);
    
    
    %--------------------------------------
    %fit the POSITION-FIXED AND PHASE-FIXED GABOR FUNCTION
    %--------------------------------------
    subplot(4,2,8);
    
    %plot the raw data
    PlotRawData(plot_x', plot_y', symbols{i}, max(px), Leye_resp, Leye_trials, Reye_resp, Reye_trials, Uncorr_resp, Uncorr_trials, null_rate);
    
    fixed_param_flags = zeros(6,1); %by default, all 6 parameters will vary
    fixed_param_values = zeros(6,1); %override these values and flags to fix a parameter
    %set the position and phase to be fixed to zero
    fixed_param_flags(6) = 1;
    fixed_param_values(6) = 0;
    fixed_param_flags(3) = 1;
    fixed_param_values(3) = 0;
    [pars{i},freq(i)] = gaborfit(means,raw,fixed_param_flags,fixed_param_values);
    
    x_interp = (px(1)): .01 : (px(length(px)));
    y_gabor = gaborfunc(x_interp, pars{i});
    y_gabor(y_gabor < 0) = 0;
    y_gauss =  pars{i}(1) + pars{i}(2)*exp(-0.5*((x_interp - pars{i}(3))/ pars{i}(4)).^2);
    y_sine =  pars{i}(1) + 0.5*pars{i}(2)*cos(2*pi*pars{i}(5)*(x_interp - pars{i}(3))+pars{i}(6) );
    
    % Now plot fitted curves; square the fitted values to counteract sqrt() above
    hold on;
    plot(x_interp, y_gabor, 'k-');
    plot(x_interp, y_gauss, 'k--');
    plot(x_interp, y_sine, 'k:');
    Gabor_SSE_position_phase_fixed(i) = gaborerr(pars{i});
    buff = sprintf('Zero Position and Phase (sse = %8.5f)', gaborerr(pars{i}) );
    title(buff);
    hold off;   
    
    y_fit = gaborfunc(px, pars{i});
    y_fit(y_fit < 0) = 0;
    %add a column of ones to yfit to make regress happy
    y_fit = [ones(length(y_fit),1) y_fit];
    [b, bint, r, rint, stats1{i}] = regress(py, y_fit);
    
    y_fit_raw = gaborfunc(plot_x', pars{i});
    y_fit_raw(y_fit_raw < 0) = 0;
    y_fit_raw = [ones(length(y_fit_raw),1) y_fit_raw];
    [b, bint, r, rint, stats2{i}] = regress(plot_y', y_fit_raw);
    
    % do a chi-square goodness of fit test
    [chi2_BothFixed, chiP_BothFixed] = chi2_test(plot_x, plot_y, 'gaborfunc', pars{i}, (length(pars{i})-1) );
    
    %print fit parameters to the screen
    subplot(4,2,7);
    axis([0 100 0 100]);    axis('off');
    xpos = -10; ypos = 100;
    font_size = 7;  bump_size = 8;
    line = sprintf('Base Rate: %8.4f', pars{i}(1));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Amplitude: %8.4f', pars{i}(2));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Gauss Ctr: %8.4f', pars{i}(3));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Gauss  SD: %8.4f', pars{i}(4));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Sine Freq: %8.4f', pars{i}(5));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Sin Phase: %8.4f', pars{i}(6));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('FT Freq: %8.4f', freq(i));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Means: Rsq=%6.3f, P=%8.6f', stats1{i}(1), stats1{i}(3));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Raw:   Rsq=%6.3f, P=%8.6f', stats2{i}(1), stats2{i}(3));
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Fit test: chi2=%6.3f, chiP=%8.6f', chi2_BothFixed, chiP_BothFixed);
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;

    names_BothFixedGabor =  sprintf('Rbase4\t Amp4\t Gctr4\t Gwid4\t Sfreq4\t Sphas4\t FTfrq4\t R2mn4\t Pmean4\t\t R2raw4\t Praw4\t\t chisq4\t chiP4\t\t');
    outstr_BothFixedGabor = sprintf('%6.2f\t %6.2f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.4f\t %10.8f\t %6.4f\t %10.8f\t %6.2f\t %10.8f\t', ...
        pars{i}, freq(i), stats1{i}(1), stats1{i}(3), stats2{i}(1), stats2{i}(3), chi2_BothFixed, chiP_BothFixed);
    
    %-------------------------------------------------------------------------------------------------------
    %do Sequential F-tests
    nfree_FullGabor = 6;
    nfree_PosFixedGabor = 5;
    nfree_PhaseFixedGabor = 5;
    nfree_BothFixedGabor = 4;
    F_Gabor_Position = ( (Gabor_SSE_pos_fixed(i) - Gabor_SSE(i))/(nfree_FullGabor-nfree_PosFixedGabor) ) / ( Gabor_SSE(i)/(length(plot_x)-nfree_FullGabor) );
    P_Gabor_Position = 1 - fcdf(F_Gabor_Position, (nfree_FullGabor-nfree_PosFixedGabor), (length(plot_x)-nfree_FullGabor) );
    F_Gabor_Phase = ( (Gabor_SSE_phase_fixed(i) - Gabor_SSE(i))/(nfree_FullGabor-nfree_PhaseFixedGabor) ) / ( Gabor_SSE(i)/(length(plot_x)-nfree_FullGabor) );
    P_Gabor_Phase = 1 - fcdf(F_Gabor_Phase, (nfree_FullGabor-nfree_PhaseFixedGabor), (length(plot_x)-nfree_FullGabor) );
    F_Gabor_Position_Phase = ( (Gabor_SSE_position_phase_fixed(i) - Gabor_SSE(i))/(nfree_FullGabor-nfree_BothFixedGabor) ) / ( Gabor_SSE(i)/(length(plot_x)-nfree_FullGabor) );
    P_Gabor_Position_Phase = 1 - fcdf(F_Gabor_Position_Phase, (nfree_FullGabor-nfree_BothFixedGabor), (length(plot_x)-nfree_FullGabor) );
    
    xpos = -20;
    ypos = ypos - bump_size;
    line = sprintf('Sequential F-tests:');
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Full vs Pos Fixed: F = %7.3f, P = %7.5f', F_Gabor_Position, P_Gabor_Position);
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Full vs Phase Fixed: F = %7.3f, P = %7.5f', F_Gabor_Phase, P_Gabor_Phase);
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Full vs Position and Phase Fixed: F = %7.3f, P = %7.5f', F_Gabor_Position_Phase, P_Gabor_Position_Phase);
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;

    names_SeqFtest = sprintf('F_Pos\t P_Pos\t\t F_Phs\t P_Phs\t\t F_Both\t P_Both\t\t ');
    outstr_SeqFtest = sprintf('%6.3f\t %10.8f\t %6.3f\t %10.8f\t %6.3f\t %10.8f\t', ...
        F_Gabor_Position, P_Gabor_Position, F_Gabor_Phase, P_Gabor_Phase, F_Gabor_Position_Phase, P_Gabor_Position_Phase);
    
end

%------------------------------------------------------------------------
%write out all relevant parameters to a cumulative text file, GCD 10/29/01
outfile = [BASE_PATH 'ProtocolSpecific\HDispTuning\GaborFit_Phase_Position.dat'];
printflag = 0;
if (exist(outfile, 'file') == 0)    %file does not yet exist
    printflag = 1;
end
fid = fopen(outfile, 'a');
if (printflag)
    fprintf(fid, '%s', [names_genl names_FullGabor names_PosFixedGabor names_PhaseFixedGabor names_BothFixedGabor names_SeqFtest]);
    fprintf(fid, '\r\n');
    printflag = 0;
end
fprintf(fid, '%s', [outstr_genl outstr_FullGabor outstr_PosFixedGabor outstr_PhaseFixedGabor outstr_BothFixedGabor outstr_SeqFtest]);
fprintf(fid, '\r\n');
fclose(fid);
%------------------------------------------------------------------------


return;


%----------------------------------------
function    PlotRawData(x, y, symb1, max_x, L_resp, L_trials, R_resp, R_trials, U_resp, U_trials, null_resp)

plot(x, y, symb1);

hold on;
errorbar(max_x*1.15, mean(L_resp), std(L_resp)/sqrt(sum(L_trials)), std(L_resp)/sqrt(sum(L_trials)), symb1);
text(max_x*1.3, mean(L_resp), 'L');
hold on;
errorbar(max_x*1.15, mean(R_resp), std(R_resp)/sqrt(sum(R_trials)), std(R_resp)/sqrt(sum(R_trials)), symb1);
text(max_x*1.3, mean(R_resp), 'R');
hold on;
errorbar(max_x*1.15, mean(U_resp), std(U_resp)/sqrt(sum(U_trials)), std(U_resp)/sqrt(sum(U_trials)), symb1);
text(max_x*1.3, mean(U_resp), 'U');


%now, get the firing rate for NULL condition trials and add spontaneous rate to plot
null_x = [min(x) max(x)];
null_y = [null_resp null_resp];
hold on;
plot(null_x, null_y, 'k--');
hold off;

return;