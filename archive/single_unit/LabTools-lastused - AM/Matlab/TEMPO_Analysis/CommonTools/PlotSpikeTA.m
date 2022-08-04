%*******************************************************************************************************************
% PlotSpikeTA.m - Generates Cross Correlograms Smoothed Shuffled Correlograms, and subtracted 
%       correlograms.  Then does bootstrapping.  In the process, subtracts shuffled from raw, and subtracts mean of 
%       subtracted correlogram to set baseline of 0.  This is done for each correlation going into the shuffles and 
%       cross correlogram
%       BJP 1/10/02

function PlotSpikeTA(data, Protocol, Analysis, SpikeChan, SpikeChan2, start_code, stop_code, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);

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
  
range = 128;		% +/- Range for generating averages
time = -range:+range; %actual time range for averages

%triangle correction factor for finite trial length
%spike_density_filter = [0.25 0.5 0.25];    
spike_density_filter = [1];

low_freq_cutoffA = 10;
high_freq_cutoffA = 100; %hz low pass cut off for integrating power

low_freq_cutoffB = 30;
high_freq_cutoffB = 70; %hz low pass cut off for integrating power

low_freq_cutoffC = 30;
high_freq_cutoffC = 200; %hz low pass cut off for integrating power

[conditions, unique_conds, param, num_conditions] = RegenerateConditionList(data, BegTrial, EndTrial, PATH, FILE, Protocol);

%now, remove trials from hor_disp and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:size(conditions,2);												% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );  

num_rows = num_conditions + 1;
num_columns = 4;
plot_offset = num_columns * 1;
figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [0 20 500 600], 'Name', 'Spike-Triggered Averages');

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
line = sprintf('Spike Density Filter: %5.3f %5.3f %5.3f %5.3f %5.3f', spike_density_filter);
text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Range For Analysis: +/-%d', range);
text(xpos,ypos, line,'FontSize',9);		ypos = ypos - bump_size;
line = sprintf('Freq Band A: %2d-%2d Hz', low_freq_cutoffA, high_freq_cutoffA );
text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Freq Band B: %2d-%2d Hz', low_freq_cutoffB, high_freq_cutoffB );
text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Freq Band C: %2d-%2d Hz', low_freq_cutoffC, high_freq_cutoffC );
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

TAve = zeros(num_conditions,2*range + 1);
Shuffled = TAve;
Subtracted = TAve;
%analyze by unique condition
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
  
    combin = 0;
    spikeTA = zeros(num_reps*num_reps, 2*range + 1);
    for trial = 1:num_reps
        start_spikebin = StartEventBin(reps(trial)) + StartOffset;
        spikes = [];
        spikes2 = [];
        spikes = data.spike_data(SpikeChan, : , reps(trial));
        spikes = find( spikes == 1); 
        spikes = spikes( (spikes > start_spikebin) & (spikes < start_spikebin + AnalInterval ) );
        
        %generate "shuffled" LFP
 
        if ~isempty(spikes)
            for trial2 = 1:num_reps
                if (neural_db2 == LFP_DB)
                    spikes2 = interp( data.lfp_data(SpikeChan2,:, reps(trial2) ), 2);                 
                else
                    spikes2 = conv(spike_density_filter, data.spike_data(SpikeChan2, :, reps(trial2))  );     
                end
                
                combin = combin + 1;    
                for spike = 1:length(spikes)  
                    spikeTA(combin,:) = spikeTA(combin,:) + spikes2(spikes(spike) - range: spikes(spike) + range);    
                end     
                spikeTA(combin,:) = spikeTA(combin,:) / length(spikes);
            end
        end
    end  
    crosses = [];
    crosses(linspace(1,(num_reps)^2, num_reps)) = 1;

    TAve(cond,:) = mean(spikeTA(logical(crosses), :) );
    Shuffled(cond,:) = mean(spikeTA(~logical(crosses),:) );
    Subtracted(cond,:) = TAve(cond,:) - Shuffled(cond, :);
    [freq, ampl] = FourierTransform_1D(0.001*(1:length(time)-1), Subtracted(cond,:), (length(time) - 1), 1,0);
    
    fourier_power(cond,:) = ampl.^2;
    power_bandA(cond) = sum(ampl((freq < high_freq_cutoffA) & (freq > low_freq_cutoffA)  ) .* ampl( (freq < high_freq_cutoffA )& (freq > low_freq_cutoffA )   )   );
    power_bandB(cond) = sum(ampl((freq < high_freq_cutoffB) & (freq > low_freq_cutoffB)  ) .* ampl( (freq < high_freq_cutoffB )& (freq > low_freq_cutoffB )   )   );
    power_bandC(cond) = sum(ampl((freq < high_freq_cutoffC) & (freq > low_freq_cutoffC)  ) .* ampl( (freq < high_freq_cutoffC )& (freq > low_freq_cutoffC )   )   );
end


RawTAmax = max(max(TAve));
RawTAmin = min(min(TAve));

STAmax = 1.5 * max(max(Subtracted));
STAmin = 1.5 * min(min(Subtracted));
Powmax = max(max(fourier_power) );

trial = 0;
for cond = 1:num_conditions
    %plot out raw spike-triggered average
    subplot(num_rows, num_columns, cond*num_columns - 3 + plot_offset);
    phandle = plot(time', Shuffled(cond,:), 'g', time', TAve(cond,:), 'r');
    xlim([-range +range]);
    ylim([RawTAmin RawTAmax]);
    if (cond == num_conditions)
        xlabel('Lag Time (ms)');    
    end
    ylabel('Raw Voltage');
    if (cond == 1)
        title('Raw STSA')    
    end
    
    %plot out spike-triggered average
    subplot(num_rows, num_columns, cond*num_columns - 2 + plot_offset);
    phandle = plot(time', Subtracted(cond,:), line_types(1));
    xlim([-range +range]);
    ylim([STAmin STAmax]);
    if (cond == num_conditions)
        xlabel('Lag Time (ms)');    
    end
    ylabel('Voltage');
    if (cond == 1)
        title('Subtracted')    
    end

    %plot out Power spectra of STAs
    subplot(num_rows, num_columns, cond*num_columns - 1 + plot_offset);
    phandle = plot(freq, fourier_power(cond,:) );
    xlim([0 200]);
    if (cond == num_conditions)
        xlabel('Frequency (Hz)');
    end
%    ylabel('Power');
    ylim([0 Powmax]);
     if (cond == 1)
        title('Power Spectra')    
    end
    
    
    %Plot identifiers and metrics
    subplot(num_rows, num_columns, cond*num_columns + plot_offset);
    axis([0 100 0 100]);                    
    axis('off');
    xpos = 0;
    ypos = 100;
    font_size = 8;
    bump_size = 30;
    line = '';
    for modality = 1: length(param) 
        line = [line param{modality}, num2str(unique_conds(cond,modality)), ' '];
    end   
    text(xpos - 30,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Power Band A: %7.5f', power_bandA(cond) );
    text(xpos,ypos, line, 'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Power Band B: %7.5f', power_bandB(cond) );
    text(xpos,ypos, line, 'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Power Band C: %7.5f', power_bandC(cond) );
    text(xpos,ypos, line, 'FontSize',font_size);		ypos = ypos - bump_size;

    axis('off');
    hold on;

end
      

