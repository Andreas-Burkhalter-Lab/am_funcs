%*******************************************************************************************************************
%	PSTH.m - This function plots PSTHs for each spike channel (in diff. colors) for the current trial condition
%			All histos are rescaled together when one of them exceeds the limits of the Y axis.
%		BJP - 4/4/03
% AnalParams = column vector containing indices for analysis
%check times for signif evoked activity bins.  Keyed off center of bins, I
%think. 4/4/03


function PSTH(data, SpikeChan, SpikeChan2, start_code, stop_code, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);

line_types = ['b-'; 'r-'; 'g-'; 'k-'; 'g-'; 'g-'; 'g-'; 'g-'];
Tempo_Defs;
ProtocolDefs;

output = 1;

%keep track of spike channels before reindexing
SpikeChanA = SpikeChan;
SpikeChanB = SpikeChan2;
num_spike_channels = size(data.spike_data,1);

switch SpikeChan
    case num_spike_channels + 1%first LFP channel
        neural_db1 = LFP_DB;
        SpikeChan = 1;
    case num_spike_channels + 2%second LFP channel
        neural_db1 = LFP_DB;
        SpikeChan = 2;
    otherwise 
        neural_db1 = SPIKE_DB;
end

%assign the proper dbase index for the neural signals and the correct indices for channel 2
%allow for the possibility of using the same channel twice -- autocorrelogram

switch SpikeChan2
    case num_spike_channels + 1%first LFP channel
        neural_db2 = LFP_DB;
        SpikeChan2 = 1;
    case num_spike_channels + 2%second LFP channel
        neural_db2 = LFP_DB;
        SpikeChan2 = 2;
    otherwise 
        neural_db2 = SPIKE_DB;
end

%need uniform duration of analysis across ALL trials - use minimum difference between start and stop 
%event bins and add offsets to find interval of analysis
AnalInterval = min( (StopEventBin + StopOffset) - (StartEventBin + StartOffset) )  ;
PostAnalInterval = min ( [mod(AnalInterval,500) (500 - mod(AnalInterval,500) )] );

%first, we'll need the bin_width (in sec) of our spike raster + events log
h = data.htb_header{neural_db1};	%for convenience
spike_bin_width = (h.skip + 1) / (h.speed_units / h.speed);

h = data.htb_header{EVENT_DB};	%for convenience
event_bin_width = (h.skip + 1) / (h.speed_units / h.speed);

[conditions, unique_conds, param, num_conditions] = RegenerateConditionList(data, BegTrial, EndTrial, PATH, FILE, Protocol);

num_modality = size(unique_conds,2) -1;

%now, remove trials from hor_disp and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:size(conditions,2);												% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );  

start_eventbin = 1;   
stop_eventbin = size(data.event_data, 2);   
%align_eventbin = find(data.event_data(1,:,1) == start_code);    

% convert start and stop times (for data range) to spike bins
total_spike_bins = AnalInterval + PostAnalInterval;
bin_width = 40;  % in ms
num_reduced_bins = floor(total_spike_bins/bin_width);

raw_hist = zeros(num_conditions,total_spike_bins);
spike_rates = zeros(1,num_conditions);
spike_rates2 = zeros(1,num_conditions);

%breaking segment of evoked activity into minisegments
minimum_segment_length = 100;

t = zeros(1, AnalInterval + PostAnalInterval);
start_e = t;
end_e = t;
    
hist_data = zeros(num_conditions, AnalInterval + PostAnalInterval);
hist_dataB = zeros(num_conditions, AnalInterval + PostAnalInterval);
   
selected_spikes1 = zeros(1, num_conditions);
selected_spikes2 = zeros(1, num_conditions);


%analyze by unique condition
for cond = 1: num_conditions
    %cond
    % These next few lines check for trials belonging to a single condition
    SetTrials(1,:) = (conditions(1,:) == unique_conds(cond,1) & select_trials); 	    
    for modality = 1: num_modality           	
        NextSetTrials = (conditions(modality,:) == unique_conds(cond,modality) & select_trials); 
        SetTrials = SetTrials & NextSetTrials;				
    end	 
    
    reps = find(SetTrials==1); 
    num_reps(cond) = length(reps);
    
    total_spikes1{cond} = zeros(1, num_reps(cond)  );
    total_spikes2{cond} = zeros(1, num_reps(cond)  );
    
    switch (neural_db1)
    case SPIKE_DB
        % initialize for first trial     
        start_bin = StartEventBin(reps(1)) + StartOffset;
        raw_spikes = data.spike_data(SpikeChan, start_bin : start_bin + AnalInterval - 1 + PostAnalInterval, reps(1)); 
        hist_data(cond,:) = hist_data(cond,:) + raw_spikes(1,:);
        [bins, counts] = SpikeBinner(raw_spikes(1,:), 1, bin_width, 0);
        hist_data3{cond}(1,:) = (counts'/(bins(2)/1000-bins(1)/1000));   
               
        total_spikes1{cond}(1) = sum(raw_spikes);
        for trial = 2: num_reps(cond)
            start_bin = StartEventBin(reps(trial)) + StartOffset;           
            raw_spikes(trial,:) = data.spike_data(SpikeChan, start_bin : start_bin + AnalInterval - 1 + PostAnalInterval, reps(trial));   
            total_spikes1{cond}(trial) = sum(raw_spikes(trial,:) );
            hist_data(cond,:) = hist_data(cond,:) +  raw_spikes(trial,:);
            [bins, counts] = SpikeBinner(raw_spikes(trial,:), 1, bin_width, 0);        
            hist_data3{cond}(trial,:) = (counts'/(bins(2)/1000-bins(1)/1000));   
        end
        
        %   spike_rates(cond);
        [bins, counts] = SpikeBinner(hist_data(cond,:) , 1, bin_width, 0);
        hist_data2(cond,:) = (counts'/(bins(2)/1000-bins(1)/1000))/num_reps(cond);   
        scale(cond) = max(hist_data2(cond,:)  );
        spike_rates(cond) = mean(hist_data2(cond,:)  );
        

    case LFP_DB
        % initialize for first trial     
        start_bin = StartEventBin(reps(1)) + StartOffset;
        hist_data = interp( data.lfp_data(SpikeChan,:, reps(1) ), 2);  
        hist_data = hist_data(  (start_bin ) : (start_bin + AnalInterval - 1 + PostAnalInterval) ); 
        hist_data2(cond,:) = resample(hist_data, 1, bin_width); 
        for trial = 2: num_reps(cond)
            start_bin = StartEventBin(reps(trial)) + StartOffset;
            hist_data = interp( data.lfp_data(SpikeChan,:, reps(1) ), 2);  
            hist_data = hist_data(  start_bin  : (start_bin + AnalInterval - 1 + PostAnalInterval) ); 
            hist_data2(cond,:) = hist_data2(cond,:) + resample(hist_data, 1, bin_width); 
        end
        hist_data2(cond,:) = hist_data2(cond,:)*(10000/bin_width)/num_reps(cond);
        scale(cond) = max(hist_data2(cond,:)  )  ;
    end
    
    switch (neural_db2)
    case SPIKE_DB
        % initialize for first trial     
        start_bin = StartEventBin(reps(1)) + StartOffset;
        raw_spikes = data.spike_data(SpikeChan2, start_bin  : start_bin + AnalInterval - 1 + PostAnalInterval, reps(1)); 
        total_spikes2{cond}(1) = sum(raw_spikes);
        hist_dataB(cond,:)  = raw_spikes(1,:);
        [bins, counts] = SpikeBinner(raw_spikes(1,:), 1, bin_width, 0);
        hist_data3B{cond}(1,:) = (counts'/(bins(2)/1000-bins(1)/1000));   

        for trial = 2: num_reps(cond)
            start_bin = StartEventBin(reps(trial)) + StartOffset;
            raw_spikes(trial,:) = data.spike_data(SpikeChan2, start_bin : start_bin + AnalInterval - 1 + PostAnalInterval, reps(trial));   
            total_spikes2{cond}(trial) = sum(raw_spikes(trial,:) );
            hist_dataB(cond,:) = hist_dataB(cond,:) +  raw_spikes(trial,:);
            [bins, counts] = SpikeBinner(raw_spikes(trial,:), 1, bin_width, 0);        
            hist_data3B{cond}(trial,:) = (counts'/(bins(2)/1000-bins(1)/1000));   
    
        end
        [bins, counts] = SpikeBinner(hist_dataB(cond,:) , 1, bin_width, 0);
        hist_data2B(cond,:) = (counts'/(bins(2)/1000-bins(1)/1000))/num_reps(cond);   
        scaleB(cond) = max(hist_data2B(cond,:)  )  ;
        spike_rates2(cond) = mean(hist_data2B(cond,:) );

    case LFP_DB
        % initialize for first trial     
        start_bin = StartEventBin(reps(1)) + StartOffset;
        hist_dataB = interp( data.lfp_data(SpikeChan2,:, reps(1) ), 2);  
        hist_dataB = hist_dataB(  start_bin  : (start_bin + AnalInterval - 1 + PostAnalInterval) ); 
        hist_data2B(cond,:) = resample(hist_dataB, 1, bin_width); 

        for trial = 2: num_reps(cond)
            start_bin = StartEventBin(reps(trial)) + StartOffset;
            hist_dataB = interp( data.lfp_data(SpikeChan2,:, reps(1) ), 2);  
            hist_dataB = hist_dataB(  (start_bin ) : (start_bin + AnalInterval - 1 + PostAnalInterval) ); 
            hist_data2B(cond,:) = hist_data2B(cond,:) + resample(hist_dataB, 1, bin_width); 
        end
        hist_data2B(cond,:) = hist_data2B(cond,:)/num_reps(cond);
        scaleB(cond) = max(hist_data2B(cond,:)  )  ;
    end
    
    time = bins';
    
    if Protocol > 200
        switch Protocol
            case FIX_1_23
                spont_condition = 1;
            case FIX_1_23_45
                spont_condition = 1;
            case FIX_VARY_HISTORY
                if cond < 6  % cond >= 4 == apertured conditions
                   spont_condition = 1;
                else
                   spont_condition = 6; %cond 4 apertured mask control
                end   
            case FIX_VARY_BACKCOLOR
                if cond < 4
                   spont_condition = 1;
                else
                   spont_condition = 4; %cond 6 apertured mask control
                end   
            otherwise
                spont_condition = 1;
        end    
    end    
    a = zeros(1, length(bins) );
    b = zeros(1, length(bins) );
    c = zeros(1, length(bins) );
    d = zeros(1, length(bins) );
    e = zeros(1, length(bins) );

    % select segments that are signicantly different and greater than first condition
    % (usually null or spontaneous condition)
        for bin = 1: size(hist_data3{spont_condition},2)  
            %use a one-tailed sign test to see if 
            [P, H] = signtest(hist_data3{spont_condition}(1:min([ (size(hist_data3{spont_condition},1) ) (size(hist_data3{cond},1))] ),bin), hist_data3{cond}(1:min([ (size(hist_data3{spont_condition},1) ) (size(hist_data3{cond},1))] ),bin), 0.025);    
            a(bin) = H;
            [P, H] = signtest(hist_data3B{spont_condition}(1:min([ (size(hist_data3{spont_condition},1) ) (size(hist_data3B{cond},1))] ),bin), hist_data3B{cond}(1:min([ (size(hist_data3{spont_condition},1) ) (size(hist_data3B{cond},1))] ),bin), 0.025);    
            b(bin) = H;
            c(bin) = median(hist_data3{spont_condition}(1:min([ (size(hist_data3{spont_condition},1) ) (size(hist_data3{cond},1))] ),bin)) < median(hist_data3{cond}(1:min([ (size(hist_data3{spont_condition},1) ) (size(hist_data3{cond},1))] ), bin )); 
            d(bin) = median(hist_data3B{spont_condition}(1:min([ (size(hist_data3{spont_condition},1) ) (size(hist_data3{cond},1))] ),bin)) < median( hist_data3B{cond}(1:min([ (size(hist_data3{spont_condition},1) ) (size(hist_data3{cond},1))] ), bin ) );
    
        end
    e = (a).*(b).*(c).*(d);
    
    %convert from bins to edges and to original sampling frequency
    t = zeros(1, AnalInterval + PostAnalInterval);
    sig_bins = [];
    sig_bins = find(e == 1);
    for i = 1: length(sig_bins)
        t( time(sig_bins(i)) - bin_width/2 + 1: time(sig_bins(i)) + bin_width/2 ) = 1;
    end

    %find starting and ending edges for significant periods of evoked
    %activity
    end_e(1:end-1) = t(1:end - 1) - t(2:end);
    start_e(2:end) = t(2:end) - t(1:end-1);
    start_times{cond} = find(start_e == 1);
    end_times{cond} = find(end_e == 1);
    start_times{cond} = start_times{cond} ( start_times{cond} <= AnalInterval );
    end_times{cond} = end_times{cond} ( end_times{cond} <= AnalInterval );
   
    %corrects for bins at the beginning or end of trial
    if ( length(start_times{cond}) > length(end_times{cond}) )
        end_times{cond}(1, end + 1) = AnalInterval; 
    end    
    if ( length(start_times{cond}) < length(end_times{cond}) )
        start_times{cond} = [ 1 start_times{cond}(1, 1: end)]
    end    

    %don't know if need this error checking to see if start and end times
    %are the same    
end


%now calculate firing rates based on selected segments
selected_FR1 = zeros(1, num_conditions);
selected_FR2 = zeros(1, num_conditions);

for cond = 1:num_conditions
   mini_segment_number = 1;
   mini_segment_start_times{cond} = [];
   mini_segment_end_times{cond}  = [];

   %check for minimum segment length
   for segment = 1: length(start_times{cond})
       if ( end_times{cond}(segment) - start_times{cond}(segment) ) >= minimum_segment_length 
           for min_segment = 1: fix(   (end_times{cond}(segment) - start_times{cond}(segment) ) / minimum_segment_length)
               mini_segment_start_times{cond}(mini_segment_number) = start_times{cond}(segment) + (min_segment - 1)*minimum_segment_length;
               mini_segment_end_times{cond}(mini_segment_number) = start_times{cond}(segment) + min_segment*minimum_segment_length;
               mini_segment_number = mini_segment_number + 1;                
           end
       end
   end
   if isempty(mini_segment_start_times{cond} )
       mini_segment_start_times{cond}(1) = 1;
       mini_segment_end_times{cond}(1) = AnalInterval;
   end
       
   selected_time = 0;
   for segment = 1:length(mini_segment_start_times{cond})
        selected_spikes1(cond) = selected_spikes1(cond) + sum(hist_data(cond, mini_segment_start_times{cond}(segment): mini_segment_end_times{cond}(segment)   ) );    
        selected_spikes2(cond) = selected_spikes2(cond) + sum(hist_dataB(cond, mini_segment_start_times{cond}(segment): mini_segment_end_times{cond}(segment)   ) );    
        selected_time = selected_time + (mini_segment_end_times{cond}(segment) - mini_segment_start_times{cond}(segment) );
   end 
   selected_FR1(cond) = selected_spikes1(cond) / (num_reps(cond) * selected_time / 1000); 
   selected_FR2(cond) = selected_spikes2(cond) / (num_reps(cond) * selected_time / 1000);
   X = [ones(length(total_spikes1{cond}), 1) total_spikes1{cond}'];
   b = regress(total_spikes2{cond}', X);
   Rsq_Spk_Cnt(cond) = b(2);
end


num_rows = num_conditions + 1;
num_columns = 2;
plot_offset = num_columns * 1;

curr_fig = figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [0 20 500 600], 'Name', 'PSTH Analysis');
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
%Next line is the source of stupid error
line = sprintf('Collected on: %s', data.htb_header{EVENT_DB}.date);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('TEMPO program: %s', data.htb_header{EVENT_DB}.pro_file);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('TEMPO config. file: %s', data.htb_header{EVENT_DB}.cfg_file);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
if (~isempty(data.eye_data))
    eh = data.htb_header{EYE_DB};
    line = sprintf('Eye data: %d channel(s) @ %6.1f Hz', eh.nchannels, eh.speed_units/(eh.speed * (eh.skip + 1)) );
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
end
if (~isempty(data.spike_data))
    sh = data.htb_header{SPIKE_DB};
    line = sprintf('Spike data: %d channel(s) @ %6.1f Hz', sh.nchannels, sh.speed_units/(sh.speed * (sh.skip + 1)) );
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
end
evh = data.htb_header{EVENT_DB};
line = sprintf('Event data: %d channel(s) @ %6.1f Hz', evh.nchannels, evh.speed_units/(evh.speed * (evh.skip + 1)) );
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Good Trials: %d out of %d total', size(data.event_data, 3), evh.sweep);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;


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
line = sprintf('Analysis Bin Width: %d ms', bin_width);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Minimum Segment Outputed: %d ms', minimum_segment_length);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;



HIST_SCALE = max([scale scaleB]);


for cond = 1: num_conditions
    subplot(num_rows, num_columns, cond*num_columns - 1 + plot_offset);
            phandle = plot(time', hist_data2(cond,:), 'k.', time', hist_data2B(cond,:), 'r.' );
            hold on;
            xx = [start_times{cond}' start_times{cond}'];
            yy = [-HIST_SCALE HIST_SCALE];
            if ~isempty(xx)
                plot(xx',yy', 'k--')
            end
            
            xx = [end_times{cond}' end_times{cond}'];
            yy = [-HIST_SCALE HIST_SCALE];
            if ~isempty(xx)
                plot(xx',yy', 'g-.')
            end
          
            
 %           phandle = plot(time, e{cond}.*f{cond}*HIST_SCALE*0.5, 'b-');
            hold on;
            %YLim([0 HIST_SCALE]);
            XLim([time(1) time(length(time))]);
            hold on  
            %set(phandle,'LineWidth', 1.5);
            xx = [0 0];
            xx2 = [AnalInterval AnalInterval];
            yy = [-HIST_SCALE HIST_SCALE];
            hold on;
            plot(xx, yy, 'b-');
            plot(xx2, yy,'b--');
            if (neural_db1 == SPIKE_DB)
                YLim([0 HIST_SCALE]);
            end
            XLim([time(1) time(length(time))]);

    hold on;
    subplot(num_rows, num_columns, cond*num_columns + plot_offset);
    axis([0 100 0 100]);                    
    axis('off');
    xpos = 0;
    ypos = 70;
    font_size = 8;
    bump_size = 30;
        
    line = '';
    for modality = 1: length(param) 
        line = [line param{modality}, num2str(unique_conds(cond,modality)), ' '];
    end   
    text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Mean FR: %4.1f spikes/s, %4.1f spikes/s', spike_rates(cond), spike_rates2(cond) );
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('FR, Selected Segments: %4.1f spikes/s, %4.1f spikes/s', selected_FR1(cond), selected_FR2(cond) );
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Number of Selected Segments: %2d', length(mini_segment_start_times{cond}) );
    text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
    

end



    %output tuning curve metrics
    if (output == 1)
        i = size(PATH,2) - 1;
        while PATH(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
            i = i - 1;
        end   
        PATHOUT = [PATH(1:i) 'Analysis\PSTH\'];
        i = size(FILE,2) - 1;
        while FILE(i) ~='.'
            i = i - 1;
        end
        FILEOUT = [FILE(1:i) 'psh'];
        eval(['save ' PATHOUT FILEOUT  ' Rsq_Spk_Cnt AnalInterval PostAnalInterval hist_data hist_dataB total_spikes1 total_spikes2 selected_spikes1 selected_spikes2 selected_FR1 selected_FR2 minimum_segment_length time hist_data2 hist_data2B spike_rates spike_rates2 mini_segment_start_times mini_segment_end_times SpikeChanA SpikeChanB neural_db1 neural_db2 start_code stop_code BegTrial EndTrial StartOffset StopOffset StartEventBin StopEventBin PATH FILE'])
        FILEOUT = [FILE(1:i) 'fig'];        
        saveas (curr_fig, [PATHOUT FILEOUT]);
    end

    
    %2 things to fix - preanal interval
    %changing from bins to edges for significant points.
    