%-----------------------------------------------------------------------------------------------------------------------
%-- SpeedHists.m -- Plots a PSTH for each speed
%--	GCD, 6/23/01
%-----------------------------------------------------------------------------------------------------------------------
function SpeedHists(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

Path_Defs;
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01
TEMPO_Defs;

%get the column of values of speeds in the dots_params matrix
speed = data.dots_params(DOTS_SPEED,:,PATCH1);

%now, get the firing rates for all the trials 
spike_rates = data.spike_rates(SpikeChan, :);

%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (speed == data.one_time_params(NULL_VALUE)) );

unique_speed = munique(speed(~null_trials)');

%now, remove trials from speed and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(speed);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

%   [speed' null_trials' select_trials']

figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [250 50 500 573], 'Name', 'Speed Histograms');
subplot(2, 1, 2);

%first, we'll need the bin_width (in sec) of our spike raster + events log
h = data.htb_header{SPIKE_DB};	%for convenience
spike_bin_width = (h.skip + 1) / (h.speed_units / h.speed);

h = data.htb_header{EVENT_DB};	%for convenience
event_bin_width = (h.skip + 1) / (h.speed_units / h.speed);

%for each different speed, compile the responses across repetitions into a PSTH and then plot it
for i=1:length(unique_speed)
    speed_select = find( (speed == unique_speed(i)) ); %trial #s
    n_reps = length(speed_select);
    
    M = [];
    for j=1:n_reps
        start_eventbin = find(data.event_data(1,:,speed_select(j)) == StartCode);
        stop_eventbin = find(data.event_data(1,:,speed_select(j)) == StopCode);
        
        % convert to spike bins, add offsets (already in spike bin format)
        start_spikebin = floor (start_eventbin*(event_bin_width/spike_bin_width)) + StartOffset;
        stop_spikebin = floor (stop_eventbin*(event_bin_width/spike_bin_width)) + StopOffset;
        
        M(j,:) = data.spike_data(1,start_spikebin:stop_spikebin, speed_select(j));     
    end

    raw_hist = sum(M,1);
    [bins{i}, counts{i}] = SpikeBinner(raw_hist',1,20,-StartOffset);
    subplot(length(unique_speed),1,i);
    bar(bins{i}, counts{i});
    
    xlim([0 1500]);
    str = sprintf('%4.1f d/s', unique_speed(i));
    xl = xlim; yl = ylim;
    text(1.01*xl(2), 0.5*yl(2), str);
    if (i == 1)
        title([PATH FILE '--Speed PSTHs']);
    end
end


%now, get the firing rate for NULL condition trials and add spontaneous rate to plot
null_rate = mean(data.spike_rates(SpikeChan, null_trials & select_trials));

output = 1;
%now output parameters about the speed tuning curve
if output == 1
    i = size(PATH,2) - 1;
    while PATH(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
        i = i - 1;
    end   
    PATHOUT = [PATH(1:i) 'Analysis\Tuning\'];
    i = size(FILE,2) - 1;
    while FILE(i) ~='.'
        i = i - 1;
    end

    FILEOUT2 = [FILE(1:i) 'spd_hists'];
    fileid = [PATHOUT FILEOUT2];
    histfid = fopen(fileid, 'w');
    %fprintf(histfid,'SpdIn\tFit\tSpeed\tAvgResp\tStdErr\tSpeed2\tSpon\n');
    
    for i=1:length(bins{3})
        for j = 1:length(unique_speed)
            fprintf(histfid, '%7.1f\t', bins{j}(i));
            fprintf(histfid,'%6.1f\t', counts{j}(i));
        end
        fprintf(histfid, '\r\n');
    end
    fclose(histfid);
    
end


return;