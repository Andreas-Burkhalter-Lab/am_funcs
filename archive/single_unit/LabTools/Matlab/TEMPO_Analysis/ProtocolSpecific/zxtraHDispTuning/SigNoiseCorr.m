%-----------------------------------------------------------------------------------------------------------------------
%-- SigNoiseCorr.m -- Calculates signal and noise correlation between two spike channels. 
%--	TU, 4/22/03
%-----------------------------------------------------------------------------------------------------------------------
function [corr_coef] = SigNoiseCorr(data, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);

TEMPO_Defs;		%needed for defines like IN_T1_WIN_CD
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01
Path_Defs;

%define the two spike channels
%SpikeChan = 3;
%SpikeChan2 = 4;

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
spike_rates1 = data.spike_rates(SpikeChan, :);
spike_rates2 = data.spike_rates(SpikeChan2, :);

%get indices of monoc. and uncorrelated controls
control_trials = logical( (hor_disp == LEYE_CONTROL) | (hor_disp == REYE_CONTROL) | (hor_disp == UNCORR_CONTROL) );

unique_hdisp = munique(hor_disp(~null_trials & ~control_trials)')

%now, remove trials from hor_disp and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(hor_disp);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

% Calculate spontaneous rates before looping through so can calculate DTI
null_resp = data.spike_rates(SpikeChan, null_trials & select_trials);
null_rate = mean(null_resp);

for i=1:length(unique_speed)	%for each different speed value, plot a separate disparity tuning curve
    speed_select = logical( (speed == unique_speed(i)) );

    plot_x = hor_disp(speed_select & ~null_trials & ~control_trials & select_trials);
    plot_y1 = spike_rates1(speed_select & ~null_trials & ~control_trials & select_trials); 
    plot_y2 = spike_rates2(speed_select & ~null_trials & ~control_trials & select_trials); 
    
    figure;
    %NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
    [px, py1, perr, pmax(i), pmin(i)] = PlotTuningCurve(plot_x', plot_y1', symbols{i}, lines{i}, 1, 1);
    [px, py2, perr, pmax(i), pmin(i)] = PlotTuningCurve(plot_x', plot_y2', symbols{i}, lines{i}, 1, 1);

    corr_coef_temp = corrcoef(py1, py2);
    sig_corr(i) = corr_coef_temp(2);
    
    for j=1:length(unique_hdisp)
        select = logical(speed_select & ~null_trials & ~control_trials & select_trials & (hor_disp == unique_hdisp(j)));
        
        z_dist1 = spike_rates1(select);
        if(sum(z_dist1) ~= 0)
            z_dist1 = (z_dist1 - mean(z_dist1))/std(z_dist1);
        end
        Z_Spikes1(select) = z_dist1;
                
        z_dist2 = spike_rates2(select);
        if(sum(z_dist2) ~= 0)
            z_dist2 = (z_dist2 - mean(z_dist2))/std(z_dist2);
        end
        Z_Spikes2(select) = z_dist2;
    end
    
    Z_Spikes1 = Z_Spikes1(speed_select & ~null_trials & ~control_trials & select_trials);
    Z_Spikes2 = Z_Spikes2(speed_select & ~null_trials & ~control_trials & select_trials);
    [Z_R,Z_P]=CORRCOEF(Z_Spikes1,Z_Spikes2);
    noise_corr(i) =  Z_R(2);
end

    %------------------------------------------------------------------------
    %write out all relevant parameters to a cumulative text file, 
    %write out one line for each stimulus speed for each neuron.
    outfile = [BASE_PATH 'ProtocolSpecific\HDispTuning\HDispTuningCorr.dat'];
    printflag = 0;
    if (exist(outfile, 'file') == 0)    %file does not yet exist
        printflag = 1;
    end
    fid = fopen(outfile, 'a');
    if (printflag)
        fprintf(fid, 'FILE\t SignalCorr\t NoiseCorr\t');
        fprintf(fid, '\r\n');
        printflag = 0;
    end
    for j = 1:length(unique_speed)
        buff = sprintf('%s\t %6.3f\t %6.3f\t', FILE, sig_corr(j), noise_corr(j));
        fprintf(fid, '%s', buff);
        fprintf(fid, '\r\n');
    end
    fclose(fid);

return;