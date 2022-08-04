%*******************************************************************************************************************
% CrossCorr - Generates Cross Correlograms Smoothed Shuffled Correlograms, and subtracted 
%       correlograms.  Then does bootstrapping.  In the process, subtracts shuffled from raw, and subtracts mean of 
%       subtracted correlogram to set baseline of 0.  This is done for each correlation going into the shuffles and 
%       cross correlogram - analyzes in frequency domain
%       BJP 1/10/01

function CrossCorr(data, Protocol, Analysis, SpikeChan, SpikeChan2, start_code, stop_code, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);

line_types = ['b-'; 'r-'; 'g-'; 'k-'; 'g-'; 'g-'; 'g-'; 'g-'];
TEMPO_Defs;
ProtocolDefs;

%keep track of spike channels before reindexing
SpikeChanA = SpikeChan;
SpikeChanB = SpikeChan2;

switch SpikeChan
    case 4 %first LFP channel
        neural_db1 = LFP_DB;
        SpikeChan = 1;
    case 5 %second LFP channel
        neural_db1 = LFP_DB;
        SpikeChan = 2;
    case 1 %first Spike Channel
        neural_db1 = SPIKE_DB;
        SpikeChan = 1;
    case 2 %second Spike Channel
        neural_db1 = SPIKE_DB;
        SpikeChan = 2;
end

%assign the proper dbase index for the neural signals and the correct indices for channel 2
%allow for the possibility of using the same channel twice -- autocorrelogram
switch SpikeChan2
    case 4 %first LFP channel
        neural_db2 = LFP_DB;
        SpikeChan2 = 1;
    case 5 %second LFP channel
        neural_db2 = LFP_DB;
        SpikeChan2 = 2;
    case 1 %first Spike Channel
        neural_db2 = SPIKE_DB;
        SpikeChan2 = 1;
    case 2 %second Spike Channel
        neural_db2 = SPIKE_DB;
        SpikeChan2 = 2;
end

%first, we'll need the bin_width (in sec) of our spike raster + events log
h = data.htb_header{neural_db1};	%for convenience
neural_bin_width = (h.skip + 1) / (h.speed_units / h.speed) * 1000;  %bin width of neural signal in ms

AnalInterval = min( (StopEventBin + StopOffset) - (StartEventBin+ StartOffset) )  ;
  
rand('state',sum(100*clock));

low_freq_cutoff = 10;
high_freq_cutoff = 100; %hz low pass cut off for integrating power

spike_density_filter = [0.25 0.5 0.25];    
%spike_density_filter = [1];
trial = 0;

time_course_window = 200;           %interval of spike trains used for correlations +/- ms
time_course_increment = 50;         %interval between time course samples (ms)
corr_window = 64;                   % +/- time lag for correlations;

%triangle correction factor for finite trial length
time = -corr_window : 1 : corr_window;
triangle2 = time_course_window/1000 - abs(time)/1000;

num_modality = 0; 
[conditions, unique_conds, param, num_conditions] = RegenerateConditionList(data, BegTrial, EndTrial, PATH, FILE, Protocol);

%now, remove trials from hor_disp and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:size(conditions,2);												% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );  

figure
num_columns = 3;
plot_offset = num_columns;
num_rows = num_conditions + plot_offset/num_columns;

set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [150 100 500 573], 'Name', 'Correlation Time Course');
subplot(10, 1, 1);
    
axis([0 100 0 100]);                    
axis('off');
xpos = -10;
ypos = 180;
font_size = 8;
bump_size = 30;

temp = strcat(PATH, FILE);
temp(temp == '\') = '/';
line = sprintf('Filename: %s %s', temp);
text(xpos,ypos, line,'FontSize',9);		ypos = ypos - bump_size;
%line = sprintf('Number of Repetitions: %3d', (size(data.spike_data,3)/num_conditions) );
%text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Spike Density Filter: %5.3f %5.3f %5.3f %5.3f %5.3f', spike_density_filter);
text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Correlation Windows: +/-%3d ms', corr_window);
text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Segments of Analysis: %3d ms', time_course_window);
text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Increments of Analysis: %2d ms', time_course_increment);
text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Freq Band For Peak Power: %2d-%2d Hz', low_freq_cutoff, high_freq_cutoff );
text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;

xpos = 40;
ypos = 150;
line = sprintf('Protocol: %d (%s)', Protocol, protocol_names{Protocol+1});
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Trial range for analysis: %d -> %d', BegTrial, EndTrial);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Spike channels analyzed: %d, %d', SpikeChanA, SpikeChanB);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Analysis Start Event: %s, offset %d ms', event_names{start_code}, StartOffset);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Analysis Stop Event: %s, offset %d ms', event_names{stop_code}, StopOffset);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;


for cond = 1 : num_conditions
    cond
    % These next few lines check for trials belonging to a single condition
    SetTrials(1,:) = (conditions(1,:) == unique_conds(cond,1) & select_trials); 	    
    for modality = 1: size(conditions,1) - 1         	
        NextSetTrials = (conditions(modality,:) == unique_conds(cond,modality) & select_trials); 
      	SetTrials = SetTrials & NextSetTrials;				
    end	 
     
    reps = find(SetTrials==1); 
    num_reps = length(reps);		   
          
    spikes = [];    
    spikes2 = [];
 	% Now do cross correlograms 
    for trial = 1:num_reps
        start_spikebin = StartEventBin(reps(trial)) + StartOffset;
        switch neural_db1
            case SPIKE_DB
                spikes(trial,:) = conv(spike_density_filter, data.spike_data(SpikeChan, start_spikebin: start_spikebin + AnalInterval, reps(trial))  ); 
            case LFP_DB
                spikes(trial,:) = data.lfp_data(SpikeChan, floor(start_spikebin/neural_bin_width): floor( (start_spikebin)/neural_bin_width ) + floor( (AnalInterval)/neural_bin_width )  , reps(trial)); 
        end   
        switch neural_db2
            case SPIKE_DB
                spikes2(trial,:) = conv(spike_density_filter, data.spike_data(SpikeChan2, start_spikebin: start_spikebin + AnalInterval, reps(trial))  );     
            case LFP_DB
                spikes2(trial,:) = data.lfp_data(SpikeChan2, floor(start_spikebin/neural_bin_width): floor( (start_spikebin)/neural_bin_width ) + floor( (AnalInterval)/neural_bin_width )   , reps(trial)); 
        end   
    end   

    
    crosses = [];
    crosses(linspace(1,(num_reps)^2, num_reps)) = 1;

    timecourse_matrix = zeros( ceil( (size(spikes,2)- time_course_window)  /time_course_increment )  , 2*corr_window + 1);
    smoothed_shuffle = timecourse_matrix;

    count = 0;
    for tau = 1 : time_course_increment: size(spikes,2) - time_course_window
        count = count + 1;
            combin = 0;
            corr_matrix = zeros( num_reps*(num_reps), 2*corr_window + 1);
            for trial1 = 1 : num_reps
                for trial2 = 1 : num_reps
                    combin = combin + 1;
                    segment = xcorr(spikes(trial1, tau : tau + time_course_window), spikes2(trial2, tau : tau + time_course_window ), corr_window)./triangle2;
                    %normalize by geometric mean spike rate
                    if ( (neural_db1 == SPIKE_DB) | (neural_db2 == SPIKE_DB) )
                        norm_factor = ( (mean( spikes(trial1, tau : tau + time_course_window) ) *1000  ) * (mean( spikes2(trial2, tau : tau + time_course_window ) )*1000 ) )^0.5;
                        if (norm_factor ~= 0)
                            segment = segment / norm_factor;
                        end    
                    end    
                    corr_matrix(combin,:) = corr_matrix(combin,:) + segment;
                end
            end     
        smoothed_shuffle(count,:) = mean(corr_matrix(~logical(crosses) , :)  );

        %subtract smoothed shuffle
        corr_matrix = corr_matrix - ones(size(corr_matrix,1),1) * smoothed_shuffle(count,:);
        %subtract mean from single trial cross correlograms
        timecourse_matrix(count,:) = mean(corr_matrix(logical(crosses), :) - mean( corr_matrix(logical(crosses), :), 2 )*ones(1, size(corr_matrix,2) ) );        
        
    end   
    

    for count = 1:size(timecourse_matrix,1)
        [freq, ampl2] = FourierTransform_1D(0.001*[1 : (size(timecourse_matrix,2) - 1)], timecourse_matrix(count, 2:end), size(timecourse_matrix,2) - 1, 1, 0);
        seg_corr_peak_power(count) = sum(ampl2((freq < high_freq_cutoff) & (freq > low_freq_cutoff)  ) .* ampl2( (freq < high_freq_cutoff)& (freq > low_freq_cutoff)   )   );
        power(count,:) = ampl2.^2; 
    end
    
    subplot(num_rows, num_columns, cond*num_columns - 1 + plot_offset);
    contour( ([1:count]*time_course_increment + 0.5*time_course_window), freq, power');    
    ylim([0 2*high_freq_cutoff]);
    xlim([(time_course_increment + 0.5*time_course_window) count*time_course_increment + 0.5*time_course_window]);
    hold on;
    ylabel('Frequency');
    if (cond == 1)
        title('Power Spectra');
    end
    if (cond == num_conditions) 
        xlabel('Time (ms)');
    end
    
    subplot(num_rows, num_columns, cond*num_columns - 2 + plot_offset);
    contour(([1:count]*time_course_increment + 0.5*time_course_window), time, timecourse_matrix');
    ylim([-50 50]);
    xlim([(time_course_increment + 0.5*time_course_window) count*time_course_increment + 0.5*time_course_window]);
    hold on;
    ylabel('Lag (ms)');
    if (cond == 1)
        title ('CCG Timecourse');
    end
    if (cond == num_conditions) 
        xlabel('Time (ms)');
    end

    subplot(num_rows, num_columns, cond*num_columns + plot_offset);
    axis([0 100 0 100]);                    
    plot( ([1:count]*time_course_increment + 0.5*time_course_window),seg_corr_peak_power);
    xlim([(time_course_increment + 0.5*time_course_window) count*time_course_increment + 0.5*time_course_window]);
    ylim([0 max(seg_corr_peak_power)]);
    ylabel ('Power');        
    if (cond == 1)
        title ('Integrated Power');    
    end    
    if (cond == num_conditions) 
        xlabel('Time (ms)');
    end    

    font_size = 8;
    bump_size = 30;
    hold on;    
    line = '';
    for modality = 1: length(param) 
        line = [line param{modality}, num2str(unique_conds(cond,modality)), ' '];
    end   
    text(0,max(seg_corr_peak_power), line,'FontSize',font_size);		ypos = ypos - bump_size;
    
end