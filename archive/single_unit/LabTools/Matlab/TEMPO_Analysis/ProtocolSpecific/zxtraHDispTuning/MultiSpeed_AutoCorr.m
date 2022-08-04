%-----------------------------------------------------------------------------------------------------------------------
%-- MultiSpeed_AutoCorr.m -- Compute and analyze autocorrelation functions from disparity tuning runs done
%--     at one or more speeds.
%--	GCD, 6/13/01
%-----------------------------------------------------------------------------------------------------------------------
function MultiSpeed_AutoCorr(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01
Path_Defs;

% not implemented yet to select output
output = 0;

symbols = {'bo' 'r*' 'go' 'mo' 'b*' 'r*' 'g*' 'c*'};
lines = {'b-' 'r-' 'g-' 'm-' 'b--' 'r--' 'g--' 'c--'};

%get the column of values of horiz. disparity in the dots_params matrix
hor_disp = data.dots_params(DOTS_HDISP,:,PATCH1);

%get indices of any NULL conditions (for measuring spontaneous activity)
null_trials = logical( (hor_disp == data.one_time_params(NULL_VALUE)) );

%get the column of speed values
speed = data.dots_params(DOTS_SPEED,:,PATCH1);
unique_speed = munique(speed(~null_trials)');

%now, get the firing rates for all the trials 
spike_rates = data.spike_rates(SpikeChan, :);

%Compute the autocorrelation for all trials
max_lag = 100;
[autocorr_matrix, lags] = ComputeAutoCorrs(data, length(hor_disp), StartCode, StopCode, StartOffset, StopOffset, max_lag);

%get indices of monoc. and uncorrelated controls
control_trials = logical( (hor_disp == LEYE_CONTROL) | (hor_disp == REYE_CONTROL) | (hor_disp == UNCORR_CONTROL) );

%now, remove trials from hor_disp and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(hor_disp);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

% Calculate spontaneous rates before looping through so can calculate DTI
null_rate = mean(data.spike_rates(SpikeChan, null_trials & select_trials));

figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'AutoCorrelation of HDisp data');

%**********************************************************************************************
%**********************************************************************************************
for i=1:length(unique_speed)	
    speed_select = logical( (speed == unique_speed(i)) );
    
    %pull out, sum, and plot autocorrelation for all trials at this speed
    select_acorr = autocorr_matrix(speed_select, :, SpikeChan);
    acorr{i} = sum(select_acorr);
    acorr2plot = sum(select_acorr)/max(sum(select_acorr));
    subplot(2, 1, 1);
    hold on;
    plot(lags, acorr{i} , lines{i});
    hold off;
    
end    
titl = sprintf('%s Autocorrelation', [PATH FILE]);
titl(titl == '\') = '/';
title(titl);
xlabel('Time Lag (ms)');
ylabel('# Coincidences');

if (unique_speed(1) == 0)
    %Take the 1D fourier transform of the autocorrelogram for the stationary data.
    %This yields the power spectrum.
    [freq, ampl] = FourierTransform_1D((lags/1000), acorr{1}, length(lags), 1, 0);
    
    subplot(2, 1, 2);
    plot(freq, ampl, 'k-');
    
    M_acorr = [lags' acorr{i}'];
    M_spec = [freq' ampl']
    
    f100 = find((freq > 97) & (freq < 103));
    ampl_100 = mean( ampl(f100) );

    %f300_500 = find(freq > 300);
    %ampl_norm = mean( ampl(f300_500) ); 

    flow = find((freq > 50) & (freq < 80));
    fhigh = find((freq > 120) & (freq < 150));
    ampl_norm = 0.5*(mean(ampl(flow) + mean(ampl(fhigh)))); 
    
    titl2 = sprintf('Power Spectrum, Static dots (f100/fnorm = %8.5f)', ampl_100/ampl_norm);
    title(titl2);
    xlabel('Frequency (Hz)');
    ylabel('Power (DC subtracted before FT)');
    
    line = sprintf('%s %12.3f %12.3f %7.5f', FILE, ampl_100, ampl_norm, ampl_100/ampl_norm);
    
    outfile = [BASE_PATH 'ProtocolSpecific\HDispTuning\HDisp_AutoCorr.dat'];
    printflag = 0;
    if (exist(outfile, 'file') == 0)    %file does not yet exist
        printflag = 1;
    end
    fid = fopen(outfile, 'a');
    if (printflag)
        fprintf(fid, 'FILE          f100    fnorm    ratio');
        fprintf(fid, '\r\n');
    end
    fprintf(fid, '%s', [line]);
    fprintf(fid, '\r\n');
    fclose(fid);
    
end

%now, print out some useful information in the upper subplot
%PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

return;