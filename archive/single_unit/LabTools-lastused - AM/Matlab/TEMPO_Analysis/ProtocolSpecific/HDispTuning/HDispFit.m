%-----------------------------------------------------------------------------------------------------------------------
%-- HDispFit.m -- Fits a disparity tuning curve, with a Gabor, sine-wave, or Gaussian function, depending
%-- on the value passed to fit_type.  NOTE: responses are square-rooted in the error functions, to help homogenize
%-- variance (a la Cumming).  A constrained version of Levenberg-Marquardt optimization is used, as
%-- implemented using 'fmincon'.
%--	GCD, 1/26/00
%-----------------------------------------------------------------------------------------------------------------------
function HDispFit(fit_type, data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01
Path_Defs;
GABOR=1;
SINE=2;
GAUSSIAN=3;
ALL=4;

% not implemented yet to select output
output = 0;

%if this flag is set, the uncorrelated responses will be included in the disparity fits.  Added by GCD, 8/13/02
USE_UNCORR_IN_FIT = 0;

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

%get the RF size for determining the disparities at which to place the Uncorr response values
RFDiam = data.one_time_params(RF_DIAMETER);

figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'HDisp Tuning Curve w/ Best Fit');

for i=length(unique_speed):length(unique_speed)	%for each different speed value, plot a separate disparity tuning curve
    %NOTE; right now this is kludged to analyze only the largest (e.g. non-zero) speed.  It may do crazy
    %things otherwise in terms of plotting.  GCD, 1/31/01
    speed_select = logical( (speed == unique_speed(i)) );
    
    plot_x = hor_disp(speed_select & ~null_trials & ~control_trials & select_trials);
    
    %TEMPORARILY convert disparity to depth, GCD 7/11/02
    %plot_x = (3.5*57./(3.5 - plot_x) - 57);  %disp in cm to depth in cm
    %plot_x = plot_x.*(1.6/48);  %scale back into range of disps so initial conditions are not terrible.
    
    plot_y = spike_rates(speed_select & ~null_trials & ~control_trials & select_trials); 
    
    %compute the Disparity Discrimination Index
    [DDI(i), var_term] = Compute_DDI(plot_x, plot_y);
    
    %compute the tuning curve Asymmetry Index
    ASI(i) = Compute_ASI(plot_x, plot_y);
    
    %NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
    [px, py, perr, pmax(i), pmin(i)] = PlotTuningCurve(plot_x', plot_y', symbols{i}, lines{i}, 1, 0);
    
    %Compute DTI from spline fit
    DTI(i) = 1 - (pmin(i).y - null_rate)/(pmax(i).y - null_rate);
    
    %store mean rates for output
    mean_rates(i, : ) = py';
    
    p_value(i) = spk_anova(plot_y, plot_x, px);
    avg_resp(i) = mean(plot_y);      
    
    
    %now, compute the monoc. and uncorrelated control values
    Leye_trials = logical( (hor_disp == LEYE_CONTROL) );
    Leye_resp = spike_rates(speed_select & select_trials & Leye_trials);
    Reye_trials = logical( (hor_disp == REYE_CONTROL) );
    Reye_resp = spike_rates(speed_select & select_trials & Reye_trials);
    Uncorr_trials = logical( (hor_disp == UNCORR_CONTROL) );
    Uncorr_resp = spike_rates(speed_select & select_trials & Uncorr_trials);
    cont_rate(i).left = mean(Leye_resp);
    cont_rate(i).right = mean(Reye_resp);
    cont_rate(i).uncorr = mean(Uncorr_resp);

    %compute the disparities to assign to the Uncorrelated responses, to include them in the fitting.
    %we set the disparity for these equal to the aperture size, since disparities larger than this would have to
    %be the same as binocularly uncorrelated dots.
    Uncorr_disp_near = -ones(1,length(Uncorr_resp))*RFDiam;
    Uncorr_disp_far = ones(1,length(Uncorr_resp))*RFDiam;
    
    means = [px py];
    raw = [plot_x' plot_y'];
    
    fit_x_uncorr = [Uncorr_disp_near plot_x Uncorr_disp_far];
    fit_y_uncorr = [Uncorr_resp plot_y Uncorr_resp];
    raw_with_uncorr = [fit_x_uncorr' fit_y_uncorr'];
    
    % do the ANOVA on sqrt(firing rate) a la Prince et al.
    [ANOVAP(i), MSgroup(i), MSerror(i)] = spk_anova_F(sqrt(plot_y), plot_x, px);
    
    fit_type = ALL;
    if ( (fit_type == GABOR) | (fit_type == ALL) )  % Fit data to gabor function
        
        if (fit_type == ALL)
            subplot(3,2,2);
        else
            subplot(2, 1, 2);
        end
        %plot the raw data
        PlotRawData(plot_x', plot_y', symbols{i}, max(px), Leye_resp, Leye_trials, Reye_resp, Reye_trials, Uncorr_resp, Uncorr_trials, null_rate);
        
        fixed_param_flags = zeros(6,1); %by default, all 6 parameters will vary
        fixed_param_values = zeros(6,1); %override these values and flags to fix a parameter
        if (USE_UNCORR_IN_FIT)
            [pars{i},freq(i)] = gaborfit(means,raw_with_uncorr,fixed_param_flags,fixed_param_values);
        else
            [pars{i},freq(i)] = gaborfit(means,raw,fixed_param_flags,fixed_param_values);
        end
        
        x_interp = (px(1)): .002 : (px(length(px)));
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
        buff = sprintf('Gabor fit (sse = %8.5f)', gaborerr(pars{i}) );
        title(buff);
        hold off;   
        
        %find the peak, trough, and 'maximum interaction position' from the Gabor fit.
        [Gabor_peak_resp, peak_i] = max(y_gabor);
        Gabor_peak_disp = x_interp(peak_i);
        [Gabor_trough_resp, trough_i] = min(y_gabor);
        Gabor_trough_disp = x_interp(trough_i);
        if (abs(Gabor_peak_resp-pars{i}(1))>abs(Gabor_trough_resp-pars{i}(1)))
            GaborMaxIntPos = Gabor_peak_disp;
        else
            GaborMaxIntPos = Gabor_trough_disp;
        end
        
        y_fit = gaborfunc(px, pars{i});
        y_fit(y_fit < 0) = 0;
        %add a column of ones to yfit to make regress happy
        y_fit = [ones(length(y_fit),1) y_fit];
        [b, bint, r, rint, stats1{i}] = regress(py, y_fit);
        
        y_fit_raw = gaborfunc(plot_x', pars{i});
        y_fit_raw(y_fit_raw < 0) = 0;
        y_fit_raw = [ones(length(y_fit_raw),1) y_fit_raw];
        [b, bint, r, rint, stats2{i}] = regress(plot_y', y_fit_raw);

        if (USE_UNCORR_IN_FIT)
            [chi2, chiP] = Chi2_Test(fit_x_uncorr, fit_y_uncorr, 'gaborfunc', pars{i}, (length(pars{i})-1) );
        else
            [chi2, chiP] = Chi2_Test(plot_x, plot_y, 'gaborfunc', pars{i}, (length(pars{i})-1) );
        end

        null_x = [min(px) max(px)];
        null_y = [null_rate null_rate];

        %----------------------------------------------------------------------------
        %also write out data in form suitable for plotting tuning curve with Origin.
        qq = size(PATH,2) - 1;
        while PATH(qq) ~='\'	%Analysis directory is one branch below Raw Data Dir
            qq = qq - 1;
        end   
        PATHOUT = [PATH(1:qq) 'Analysis\Tuning\'];
        qq = size(FILE,2) - 1;
        while FILE(qq) ~='.'
            qq = qq - 1;
        end
        FILEOUT2 = [FILE(1:qq) 'hdsp_Gabor_fit'];
        fileid = [PATHOUT FILEOUT2];
        proffid = eval(['fopen(fileid, ''w'')']);
    
        fprintf(proffid,'HDspI\tGbFit\tHDisp\tAvgRsp\tStdErr\tMonDsp\tMonoc\tMonErr\tLabel\tSpDsp\tSpon\n');
        for mm=1:length(x_interp)
            fprintf(proffid,'%6.3f\t%6.3f\t', x_interp(mm), y_gabor(mm));
            if (mm <= length(px))
                fprintf(proffid,'%6.2f\t%6.2f\t%6.3f\t', px(mm), py(mm), perr(mm));
            else
                fprintf(proffid,'\t\t\t');
            end
            if (mm == 1)
                fprintf(proffid,'%6.2f\t%6.2f\t%6.2f\t%6.1f\t%6.2f\t%6.2f\n',1.15*px(length(px)),mean(Leye_resp),std(Leye_resp)/sqrt(sum(Leye_trials)),-99,null_x(mm),null_y(mm));
            elseif (mm==2)		
                fprintf(proffid,'%6.2f\t%6.2f\t%6.2f\t%6.1f\t%6.2f\t%6.2f\n',1.15*px(length(px)),mean(Reye_resp),std(Reye_resp)/sqrt(sum(Reye_trials)),99,null_x(mm),null_y(mm));
            elseif (mm==3)		
                fprintf(proffid,'%6.2f\t%6.2f\t%6.2f\t%6.1f\t%6.2f\t%6.2f\n',1.15*px(length(px)),mean(Uncorr_resp),std(Uncorr_resp)/sqrt(sum(Uncorr_trials)),98,null_x(1),null_y(1));
            else
                fprintf(proffid,'\t\t\t\t\t\n');
            end
        end
    
        fclose(proffid);
        %----------------------------------------------------------------------------
        
        %print fit parameters to the screen
        if (fit_type == ALL)
            subplot(3,2,1);
            axis([0 100 0 100]);    axis('off');
            xpos = -10; ypos = 110;
            font_size = 8;  bump_size = 9;
            temp = strcat(PATH, FILE);
            temp(temp == '\') = '/';
            % this prevents a stupid error from appearing on the screen
            line = sprintf('File: %s', temp);
            text(xpos-15,ypos, line,'FontSize',font_size+1);		ypos = ypos - bump_size;
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
            line = sprintf('ANOVA P:   %10.6f', ANOVAP(i));
            text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
        end

        % put data into strings for output to file
        names_genl =  sprintf('FILE\t\t PrDir\t PrSpd\t PrHDsp\t RFx\t RFy\t RFdiam\t Lrsp\t Lse\t Rrsp\t Rse\t Ursp\t Use\t GbPkR\t GbPkD\t GbTrR\t GbTrD\t MxIntP\t AnovP\t\t MSgrp\t\t MSerr\t\t');
        outstr_genl = sprintf('%s\t %6.1f\t %6.2f\t %6.3f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.3f\t %6.2f\t %6.3f\t %6.3f\t %10.8f\t %10.3f\t %10.3f\t', ...
            FILE, data.neuron_params(PREFERRED_DIRECTION, 1), data.neuron_params(PREFERRED_SPEED, 1), data.neuron_params(PREFERRED_HDISP, 1), data.neuron_params(RF_XCTR, 1), ... 
            data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1), mean(Leye_resp),std(Leye_resp)/sqrt(sum(Leye_trials)), ...
            mean(Reye_resp),std(Reye_resp)/sqrt(sum(Reye_trials)), mean(Uncorr_resp),std(Uncorr_resp)/sqrt(sum(Uncorr_trials)), ... 
            Gabor_peak_resp, Gabor_peak_disp, Gabor_trough_resp, Gabor_trough_disp, GaborMaxIntPos, ANOVAP(i), MSgroup(i), MSerror(i));
        
        [phase360, phase180, phase90] = AngleWrap(pars{i}(6)*180/pi);
        if (phase360 > 180)
            phase360 = phase360 - 360;  %fold into -180->180 range
        end
    
        names_Gabor =  sprintf('Rbase1\t Amp1\t Gctr1\t Gwid1\t Sfreq1\t Sphs1\t GbP360\t GbP180\t GbP90\t FTfrq1\t R2mn1\t Pmean1\t\t R2raw1\t Praw1\t\t chisq1\t chiP1\t\t');
        outstr_Gabor = sprintf('%6.2f\t %6.2f\t %6.3f\t %6.3f\t %6.3f\t %6.3f\t %6.1f\t %6.1f\t %6.1f\t %6.3f\t %6.4f\t %10.8f\t %6.4f\t %10.8f\t %6.2f\t %10.8f\t', ...
            pars{i}, phase360, phase180, phase90, freq(i), stats1{i}(1), stats1{i}(3), stats2{i}(1), stats2{i}(3), chi2, chiP);
        
    end  %(fit_type == GABOR)
    if ( (fit_type == SINE) | (fit_type == ALL) )  % Fit sinusoid to data
        
        if (fit_type == ALL)
            subplot(3,2,4);
        else
            subplot(2, 1, 2);
        end
        %plot the raw data
        PlotRawData(plot_x', plot_y', symbols{i}, max(px), Leye_resp, Leye_trials, Reye_resp, Reye_trials, Uncorr_resp, Uncorr_trials, null_rate);
        
        if (USE_UNCORR_IN_FIT)
            [pars{i}] = sine_fit(means,raw_with_uncorr);
        else
            [pars{i}] = sine_fit(means,raw);
        end
        
        x_interp = (px(1)): .01 : (px(length(px)));
        y_sine =  sine_func(x_interp,pars{i});
        y_sine(y_sine < 0) = 0;
        
        % Now plot fitted curve
        hold on;
        plot(x_interp, y_sine, 'k-');                                                                    
        Sine_SSE(i) = sine_err(pars{i});
        buff = sprintf('Sine fit (sse = %8.5f)', sine_err(pars{i}) );
        title(buff);
        hold off;   
        
        y_fit = sine_func(px, pars{i});
        y_fit(y_fit < 0) = 0;
        %add a column of ones to yfit to make regress happy
        y_fit = [ones(length(y_fit),1) y_fit];
        [b, bint, r, rint, stats1{i}] = regress(py, y_fit);
        
        y_fit_raw = sine_func(plot_x', pars{i});
        y_fit_raw(y_fit_raw < 0) = 0;
        y_fit_raw = [ones(length(y_fit_raw),1) y_fit_raw];
        [b, bint, r, rint, stats2{i}] = regress(plot_y', y_fit_raw);
        
        if (USE_UNCORR_IN_FIT)
            [chi2, chiP] = Chi2_Test(fit_x_uncorr, fit_y_uncorr, 'sine_func', pars{i}, length(pars{i}) );
        else
            [chi2, chiP] = Chi2_Test(plot_x, plot_y, 'sine_func', pars{i}, length(pars{i}) );
        end
        
        %print fit parameters to the screen
        if (fit_type == ALL)
            subplot(3,2,3);
            axis([0 100 0 100]);    axis('off');
            xpos = -10; ypos = 110;
            font_size = 9;  bump_size = 9;
            line = sprintf('Base Rate: %8.4f', pars{i}(1));
            text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
            line = sprintf('Amplitude: %8.4f', pars{i}(2));
            text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
            line = sprintf('Sine Freq: %8.4f', pars{i}(3));
            text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
            line = sprintf('Sin Phase: %8.4f', pars{i}(4));
            text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
            line = sprintf('Means: Rsq=%6.3f, P=%8.6f', stats1{i}(1), stats1{i}(3));
            text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
            line = sprintf('Raw:   Rsq=%6.3f, P=%8.6f', stats2{i}(1), stats2{i}(3));
            text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
        end

        [phase360, phase180, phase90] = AngleWrap(pars{i}(4)*180/pi);
        if (phase360 > 180)
            phase360 = phase360 - 360;  %fold into -180->180 range
        end
    
        names_Sine =  sprintf('Rbase2\t Amp2\t Sfreq2\t Sphs2\t SnP360\t SnP180\t SnP90\t R2mn2\t Pmean2\t\t R2raw2\t Praw2\t\t chisq2\t chiP2\t\t');
        outstr_Sine = sprintf('%6.2f\t %6.2f\t %6.3f\t %6.3f\t %6.1f\t %6.1f\t %6.1f\t %6.4f\t %10.8f\t %6.4f\t %10.8f\t %6.2f\t %10.8f\t', ...
            pars{i}, phase360, phase180, phase90, stats1{i}(1), stats1{i}(3), stats2{i}(1), stats2{i}(3), chi2, chiP);
    end
    if ( (fit_type == GAUSSIAN) | (fit_type == ALL) )  % Fit Gaussian to data
        
        if (fit_type == ALL)
            subplot(3,2,6);
        else
            subplot(2, 1, 2);
        end
        %plot the raw data
        PlotRawData(plot_x', plot_y', symbols{i}, max(px), Leye_resp, Leye_trials, Reye_resp, Reye_trials, Uncorr_resp, Uncorr_trials, null_rate);
        
        if (USE_UNCORR_IN_FIT)
            [pars{i}] = gaussfit(means,raw_with_uncorr,1);
        else
            [pars{i}] = gaussfit(means,raw,1);
        end
        
        x_interp = (px(1)): .01 : (px(length(px)));
        y_gauss =  gaussfunc(x_interp,pars{i});
        y_gauss(y_gauss < 0) = 0;
        
        % Now plot fitted curve
        hold on;
        plot(x_interp, y_gauss, 'k-');                                                                    
        Gaussian_SSE(i) = gausserr(pars{i});
        buff = sprintf('Gaussian fit (sse = %8.5f)', gausserr(pars{i}) );
        title(buff);
        hold off;   
        
        y_fit = gaussfunc(px, pars{i});
        y_fit(y_fit < 0) = 0;
        %add a column of ones to yfit to make regress happy
        y_fit = [ones(length(y_fit),1) y_fit];
        [b, bint, r, rint, stats1{i}] = regress(py, y_fit);
        
        y_fit_raw = gaussfunc(plot_x', pars{i});
        y_fit_raw(y_fit_raw < 0) = 0;
        y_fit_raw = [ones(length(y_fit_raw),1) y_fit_raw];
        [b, bint, r, rint, stats2{i}] = regress(plot_y', y_fit_raw);

        if (USE_UNCORR_IN_FIT)
            [chi2, chiP] = Chi2_Test(fit_x_uncorr, fit_y_uncorr, 'gaussfunc', pars{i}, length(pars{i}) );
        else
            [chi2, chiP] = Chi2_Test(plot_x, plot_y, 'gaussfunc', pars{i}, length(pars{i}) );
        end
        
        %print fit parameters to the screen
        if (fit_type == ALL)
            subplot(3,2,5);
            axis([0 100 0 100]);    axis('off');
            xpos = -10; ypos = 110;
            font_size = 9;  bump_size = 9;
            line = sprintf('Base Rate: %8.4f', pars{i}(1));
            text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
            line = sprintf('Amplitude: %8.4f', pars{i}(2));
            text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
            line = sprintf('Gauss Ctr: %8.4f', pars{i}(3));
            text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
            line = sprintf('Gauss  SD: %8.4f', pars{i}(4));
            text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
            line = sprintf('Means: Rsq=%6.3f, P=%8.6f', stats1{i}(1), stats1{i}(3));
            text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
            line = sprintf('Raw:   Rsq=%6.3f, P=%8.6f', stats2{i}(1), stats2{i}(3));
            text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
        end

        names_Gauss =  sprintf('Rbase3\t Amp3\t Gctr3\t Gwid3\t R2mn3\t Pmean3\t\t R2raw3\t Praw3\t\t chisq3\t chiP3\t\t');
        outstr_Gauss = sprintf('%6.2f\t %6.2f\t %6.3f\t %6.3f\t %6.4f\t %10.8f\t %6.4f\t %10.8f\t %6.2f\t %10.8f\t', ...
            pars{i}, stats1{i}(1), stats1{i}(3), stats2{i}(1), stats2{i}(3), chi2, chiP);
        
    end

    %do Sequential F-tests
    if (fit_type == ALL)
        nfree_Gabor = 6;
        nfree_Gauss = 4;
        nfree_Sine = 4;
        F_Gabor_Sine = ( (Sine_SSE(i) - Gabor_SSE(i))/(nfree_Gabor-nfree_Sine) ) / ( Gabor_SSE(i)/(length(plot_x)-nfree_Sine) );
        P_Gabor_Sine = 1 - fcdf(F_Gabor_Sine, (nfree_Gabor-nfree_Sine), (length(plot_x)-nfree_Gabor) );
        F_Gabor_Gauss = ( (Gaussian_SSE(i) - Gabor_SSE(i))/(nfree_Gabor-nfree_Gauss) ) / ( Gabor_SSE(i)/(length(plot_x)-nfree_Gabor) );
        P_Gabor_Gauss = 1 - fcdf(F_Gabor_Gauss, (nfree_Gabor-nfree_Gauss), (length(plot_x)-nfree_Gabor) );
        
        ypos = ypos - bump_size;
        line = sprintf('Sequential F-tests:');
        text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
        line = sprintf('Gabor - Sine: F = %7.3f, P = %7.5f', F_Gabor_Sine, P_Gabor_Sine);
        text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
        line = sprintf('Gabor - Gauss: F = %7.3f, P = %7.5f', F_Gabor_Gauss, P_Gabor_Gauss);
        text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
        
        names_SeqFtest = sprintf('F_GbSn\t P_GbSn\t\t F_GbGs\t P_GbGs\t\t ');
        outstr_SeqFtest = sprintf('%6.3f\t %10.8f\t %6.3f\t %10.8f\t', F_Gabor_Sine, P_Gabor_Sine, F_Gabor_Gauss, P_Gabor_Gauss);
        
    end

end

yl = YLim;
YLim([0 yl(2)]);	% set the lower limit of the Y axis to zero
XLabel('Horizontal Disparity(deg)');
YLabel('Response (spikes/sec)');

%now, print out some useful information in the upper subplot
if (fit_type ~= ALL)
    subplot(2, 1, 1);
    PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);
end

%now print out useful values for H.Disp specific 
% pmax, pmin, py, 
if (fit_type ~= ALL)
    PrintHDispData(p_value, avg_resp, pmax, pmin, px, null_rate, cont_rate, unique_speed, PATH, FILE, DDI, DTI, corr_coef, ASI);
end

%------------------------------------------------------------------------
%write out all relevant parameters to a cumulative text file, GCD 10/29/01
outfile = [BASE_PATH 'ProtocolSpecific\HDispTuning\HDispFit_GaborSineGauss.dat'];
printflag = 0;
if (exist(outfile, 'file') == 0)    %file does not yet exist
    printflag = 1;
end
fid = fopen(outfile, 'a');
if (printflag)
    fprintf(fid, '%s', [names_genl names_Gabor names_Sine names_Gauss names_SeqFtest]);
    fprintf(fid, '\r\n');
    printflag = 0;
end
fprintf(fid, '%s', [outstr_genl outstr_Gabor outstr_Sine outstr_Gauss outstr_SeqFtest]);
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