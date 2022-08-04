%Edited on 1/20/08 by RLS to incorporate a seperate 'Pooling' function for use
%with calculating latency
%Edited on 10/05/07 by RLS to display file name on the top of the plots
%Edited on 06/19/06 by CMA to remove all references to SpikeChan2 or
%spikes2, which produce the error "??? Undefined function or variable
%'SpikeChan2'." when run with our data.  The original file is backed up in
%the folder Chris's Lab Tools Backup located on the desktop of this
%machine.

function SF_PSTH(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol, Pooling);
line_types = ['b-'; 'r-'; 'g-'; 'k-'; 'g-'; 'g-'; 'g-'; 'g-'];
Tempo_Defs;
ProtocolDefs;

TEMPO_Defs;		%needed for defines like IN_T1_WIN_CD
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01
Path_Defs;

symbols = {'ko' 'k*' 'go' 'mo' 'b*' 'r*' 'g*' 'c*'};
lines = {'k-' 'm-' 'g-' 'b-' 'b--' 'r--' 'g--' 'c--'};
colors = {'k.' 'b.'};
global pool_size pool_data1 pool_data2 pool_trials sample_trials bin_size window_size, 

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
count=1;
%subplot(2, 1, 2);
for j = 1:length(unique_stim_type);

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
            count = count+1;
            spikes =  find(data.spike_data(SpikeChan,(StartEventBin(trial) + StartOffset - 500):(StopEventBin(trial) + StopOffset + 500),trial) == 1);
            %plot (spikes/1000,count*ones(1,length(spikes)),colors{j});
            if j == 1;
                total_spikes1 = cat(2,total_spikes1,spikes);
            elseif j == 2;
                total_spikes2 = cat(2,total_spikes2,spikes);
            else
                'Error: More than two grating types'
            end
            hold on;
        end
    end

    trials(j) = count-1;

end

bin_matrix = [0:0.05:3]; %50 msec bins
bin_centers = [0.025:0.05:2.975];
bin_centers = [bin_centers, 3];
[hits1] = histc(total_spikes1/1000, bin_matrix);
if ~isempty(total_spikes2)
    [hits2] = histc(total_spikes2/1000, bin_matrix);
end

if Pooling == 1;
    pool_size = pool_size + 1;
    pool_trials = pool_trials + sample_trials;
    pool_data1 = cat(2, pool_data1, total_spikes1);
    if ~isempty(total_spikes2)
        pool_data2 = cat(2, pool_data2, total_spikes2);
    end
    if pool_size == 5;
        Pool_PSTH(bin_size, window_size, data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol); 
    end
    disp('Pooling complete');
    return
end

if plots_on == 1;
    figure;
    set(gcf,'PaperPosition', [.51 .51 20 27.5], 'Position', [450 50 500 573], 'Name', 'Optimized Spatial Frequency PSTH');
    hold on
    subplot(20, 1, 2:20);
    plot(bin_centers,hits1,lines{1})
    if ~isempty(total_spikes2)
        plot(bin_centers,hits2,lines{2})
    end
    title('Optimized Spatial Frequency PSTH')
    subplot(20, 1, 1);
    PrintFileName(FILE);
end

%Normalize the responses to a peak of 1 so that a few very active cells do
%not bias the results
normalized_hits1 = hits1/max(hits1);
if ~isempty(total_spikes2)
    normalized_hits2 = hits2/max(hits2);
end

if plots_on == 1;
    figure;
    set(gcf,'PaperPosition', [.51 .51 20 27.5], 'Position', [450 50 500 573], 'Name', 'Optimized Spatial Frequency PSTH (Normalized)');
    hold on
    subplot(20, 1, 2:20);
    plot(bin_centers,normalized_hits1,lines{1})
    if ~isempty(total_spikes2)
        plot(bin_centers,normalized_hits2,lines{2})
    end
    title('Optimized Spatial Frequency PSTH - Normalized')
    subplot(20, 1, 1);
    PrintFileName(FILE);
end


output = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non-normalized Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (output == 1)

    %------------------------------------------------------------------------
    %write out all relevant parameters to a cumulative text file, GCD 8/08/01
    %write out one line for each stimu_type for each neuron.
    outfile = [BASE_PATH 'ProtocolSpecific\GratingSpatialFreq\OptimizedSF_PSTH.dat'];
    printflag = 0;
    if (exist(outfile, 'file') == 0)    %file does not yet exist
        printflag = 1;
    end
    fid = fopen(outfile, 'a');
    if (printflag)
        fprintf(fid,...
            '%5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t', bin_centers)
        fprintf(fid, '\r\n');
        printflag = 0;
    end
    %for j = 1:length(unique_stim_type)
    fprintf(fid,...
        '%5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t', hits1)
    fprintf(fid, '\r\n');
    %end
    fclose(fid);
    %------------------------------------------------------------------------

end  %if (output == 1)

if ~isempty(total_spikes2)
    if (output == 1)

        %------------------------------------------------------------------------
        %write out all relevant parameters to a cumulative text file, GCD 8/08/01
        %write out one line for each stimu_type for each neuron.
        outfile = [BASE_PATH 'ProtocolSpecific\GratingSpatialFreq\OptimizedSF_PSTH2.dat'];
        printflag = 0;
        if (exist(outfile, 'file') == 0)    %file does not yet exist
            printflag = 1;
        end
        fid = fopen(outfile, 'a');
        if (printflag)
            fprintf(fid,...
                '%5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t', bin_centers)
            fprintf(fid, '\r\n');
            printflag = 0;
        end
        %for j = 1:length(unique_stim_type)
        fprintf(fid,...
            '%5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t', hits2)
        fprintf(fid, '\r\n');
        %end
        fclose(fid);
        %------------------------------------------------------------------------

    end  %if (output == 1)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalized Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (output == 1)

    %------------------------------------------------------------------------
    %write out all relevant parameters to a cumulative text file, GCD 8/08/01
    %write out one line for each stimu_type for each neuron.
    outfile = [BASE_PATH 'ProtocolSpecific\GratingSpatialFreq\OptimizedSF_PSTH_Norm.dat'];
    printflag = 0;
    if (exist(outfile, 'file') == 0)    %file does not yet exist
        printflag = 1;
    end
    fid = fopen(outfile, 'a');
    if (printflag)
        fprintf(fid,...
            '%5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t', bin_centers)
        fprintf(fid, '\r\n');
        printflag = 0;
    end
    %for j = 1:length(unique_stim_type)
    fprintf(fid,...
        '%5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t', normalized_hits1)
    fprintf(fid, '\r\n');
    %end
    fclose(fid);
    %------------------------------------------------------------------------

end  %if (output == 1)

if ~isempty(total_spikes2)
    if (output == 1)

        %------------------------------------------------------------------------
        %write out all relevant parameters to a cumulative text file, GCD 8/08/01
        %write out one line for each stimu_type for each neuron.
        outfile = [BASE_PATH 'ProtocolSpecific\GratingSpatialFreq\OptimizedSF_PSTH_Norm2.dat'];
        printflag = 0;
        if (exist(outfile, 'file') == 0)    %file does not yet exist
            printflag = 1;
        end
        fid = fopen(outfile, 'a');
        if (printflag)
            fprintf(fid,...
                '%5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t', bin_centers)
            fprintf(fid, '\r\n');
            printflag = 0;
        end
        %for j = 1:length(unique_stim_type)
        fprintf(fid,...
            '%5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t %5.3f\t', normalized_hits2)
        fprintf(fid, '\r\n');
        %end
        fclose(fid);
        %------------------------------------------------------------------------

    end  %if (output == 1)
end
return;














