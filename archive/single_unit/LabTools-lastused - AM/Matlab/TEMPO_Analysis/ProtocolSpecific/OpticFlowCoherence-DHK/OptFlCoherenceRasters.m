%Edited on 06/19/06 by CMA to remove all references to SpikeChan2 or
%spikes2, which produce the error "??? Undefined function or variable
%'SpikeChan2'." when run with our data.  The original file is backed up in
%the folder Chris's Lab Tools Backup located on the desktop of this
%machine.

function OptFlCoherenceRasters(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
line_types = ['b-'; 'r-'; 'g-'; 'k-'; 'g-'; 'g-'; 'g-'; 'g-'];
Tempo_Defs;
ProtocolDefs;

TEMPO_Defs;		%needed for defines like IN_T1_WIN_CD
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01
Path_Defs;

symbols = {'ko' 'k*' 'go' 'mo' 'b*' 'r*' 'g*' 'c*'};
lines = {'k-' 'k--' 'g-' 'm-' 'b--' 'r--' 'g--' 'c--'};
colors = {'k.' 'b.'};

%get the column of values of SFs in gratings_params matrix
coher = data.dots_params(DOTS_COHER,:,PATCH1);

%get indices of any NULL conditions (for measuring spontaneous activity)
null_trials = logical( (coher == data.one_time_params(NULL_VALUE)) );

%get the column of stimulus type values
direc = data.dots_params(DOTS_DIREC,:,PATCH1);
unique_direc = munique(direc(~null_trials)');

%now, get the firing rates for all the trials
spike_rates = data.spike_rates(SpikeChan, :);

unique_coher = munique(coher(~null_trials)');

%now, remove trials that do not fall between BegTrial and EndTrial
trials = 1:length(coher);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );


total_spikes1 = [];
total_spikes2 = [];
count=1;
figure(2);
set(gcf,'PaperPosition', [.51 .51 20 27.5], 'Position', [450 50 500 573], 'Name', 'Optimized Coherence Raster');
for j = 1:length(unique_direc);

    direc_select = logical( (direc == unique_direc(j)) );
    
    plot_x=coher(~null_trials & select_trials & direc_select);
    plot_y=spike_rates(~null_trials & select_trials & direc_select);

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
    for trial=1:length(coher);
        if (coher(trial) == optimized_stimulus & direc(trial) == unique_direc(j));
            count = count+1;
            spikes =  find(data.spike_data(SpikeChan,(StartEventBin(trial) + StartOffset - 500):(StopEventBin(trial) + StopOffset + 500),trial) == 1);
            plot (spikes/1000,count*ones(1,length(spikes)),colors{j});
            if j == 1;
                total_spikes1 = cat(2,total_spikes1,spikes);
            elseif j == 2;
                total_spikes2 = cat(2,total_spikes2,spikes);
            else
                'Error: More than two directions'
            end
            hold on;
        end
    end
    
    trials(j) = count-1;
    
end

stim_bar = [1:3000];
on_off = ones(size(stim_bar));
on_off(1:500)=0;
on_off(2500:3000)=0;
plot(stim_bar/1000,on_off,'k')
axis([0 3 0 (count+1)])
xlabel('Time (s)')
ylabel('Trial')
title('Optimized Coherence Spike Rasters')