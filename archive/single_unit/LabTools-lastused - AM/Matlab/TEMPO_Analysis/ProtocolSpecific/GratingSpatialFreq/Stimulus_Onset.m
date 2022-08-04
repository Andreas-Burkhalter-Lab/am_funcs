%Edited on 8/6/07 by CMA.  This program works using the method from Hanes
%et al 1995.  

function Stimulus_Onset(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
line_types = ['b-'; 'r-'; 'g-'; 'k-'; 'g-'; 'g-'; 'g-'; 'g-'];
Tempo_Defs;
ProtocolDefs;

TEMPO_Defs;		%needed for defines like IN_T1_WIN_CD
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01
Path_Defs;

symbols = {'ko' 'k*' 'go' 'mo' 'b*' 'r*' 'g*' 'c*'};
lines = {'k-' 'm-' 'g-' 'b-' 'b--' 'r--' 'g--' 'c--'};
colors = {'k.' 'b.'};

%get the column of values of SFs in gratings_params matrix
sf = data.gratings_params(GRAT_SPATIAL_FREQ,:,PATCH1);

%get indices of any NULL conditions (for measuring spontaneous activity)
null_trials = logical( (sf == data.one_time_params(NULL_VALUE)) );

%get the column of stimulus type values
stim_type = data.gratings_params(GRAT_TYPE,:,PATCH1);
unique_stim_type = munique(stim_type(~null_trials)');

%now, get the firing rates for all the trials
spike_rates = data.spike_rates(SpikeChan, :);

unique_sf = munique(sf(~null_trials)');

%now, remove trials that do not fall between BegTrial and EndTrial
trials = 1:length(sf);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

plots_on = 1;

total_spikes1 = [];
total_spikes2 = [];
null_resp = data.spike_rates(SpikeChan, null_trials & select_trials);
null_rate = mean(null_resp);
StartT = 501;
StopT = 3000;

for j = 1:length(unique_stim_type);
    count=1;
    grating_select = logical( (stim_type == unique_stim_type(j)) );
    
    plot_x=sf(~null_trials & select_trials & grating_select);
    plot_y=spike_rates(~null_trials & select_trials & grating_select);

    %NOTE:  Inputs to PlotTuningCurve.m must be column vectors, not row vectors
    %because of munique().
    [px, py, perr] = PlotTuningCurve(plot_x', plot_y', symbols{j}, lines{j}, 1, 0);
    %NOTE: last argument=0 means just get output, no plot
    
    %py is the response and px is the stimulus value.  We are interested in
    %the cell's response at optimized conditions, therefore, we find the
    %max of py, and select the corresponding px as stimulus value that
    %results in optimized response.
    [max_response, index] = max(py);
    optimized_stimulus = px(index);

 
    hold on
    for trial=1:length(sf);
        if (sf(trial) == optimized_stimulus & stim_type(trial) == unique_stim_type(j));
            
            %Find all the spikes between the 500ms before stimulus and the
            %500ms after stimulus
            spikes =  find(data.spike_data(SpikeChan,(StartEventBin(trial) + StartOffset - 500):(StopEventBin(trial) + StopOffset + 500),trial) == 1);
            if j == 1;
                total_spikes1 = cat(2,total_spikes1,spikes);
            elseif j == 2;
                total_spikes2 = cat(2,total_spikes2,spikes);
            else
                'Error: More than two grating types'
            end
            count = count+1;
            hold on;
        end
    end
    
    trials(j) = count-1;
    
end



total_spikes1 = sort(total_spikes1);
total_spikes2 = sort(total_spikes2);
[response_time(1),EOB(1),P(1)] = burst(total_spikes1,StartT,StopT)
if ~isempty(total_spikes2)
    [response_time(2),EOB(2),P(2)] = burst(total_spikes2,StartT,StopT)
end

output = 1;
if (output == 1)

    %------------------------------------------------------------------------
    %write out all relevant parameters to a cumulative text file, GCD 8/08/01
    %write out one line for each stimu_type for each neuron.
    outfile = [BASE_PATH 'ProtocolSpecific\GratingSpatialFreq\SF_Onset.dat'];
    printflag = 0;
    if (exist(outfile, 'file') == 0)    %file does not yet exist
        printflag = 1;
    end
    fid = fopen(outfile, 'a');
    if (printflag)
        fprintf(fid, 'FILE\t\t\t PrfSF\t GrType\t Onset\t EOB\t\t P\t');
        fprintf(fid, '\r\n');
        printflag = 0;
    end
    for j = 1:length(unique_stim_type)
        buff = sprintf('%s\t %5.2f\t %5.0f\t %5.3f\t %5.3f\t %5.3e\t', ...
            FILE, data.neuron_params(PREFERRED_SPATIAL_FREQ, 1),...
            unique_stim_type(j), response_time(j), EOB(j), P(j));
        fprintf(fid, '%s', buff);
        fprintf(fid, '\r\n');
    end
    fclose(fid);
    %------------------------------------------------------------------------
end



