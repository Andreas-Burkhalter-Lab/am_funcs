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
    case 3 %third Spike Channel
        neural_db2 = SPIKE_DB;
        SpikeChan2 = 3;
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
    case 3 %third Spike Channel
        neural_db2 = SPIKE_DB;
        SpikeChan2 = 3;
        
end

%first, we'll need the bin_width (in sec) of our spike raster + events log
h = data.htb_header{neural_db1};	%for convenience
neural_bin_width = (h.skip + 1) / (h.speed_units / h.speed) * 1000;  %bin width of neural signal in ms

AnalInterval = min( (StopEventBin + StopOffset) - (StartEventBin+ StartOffset) )  ;
  
rand('state',sum(100*clock));

num_bootstraps = 10000;

range = 128;		% Range for generating cross correlograms
disp_range = 128;    %Range for histograms
time = -range:neural_bin_width:+range; %actual time range for cross correlograms

peak_range = 5; % +/- range for width of averaging across peak
peak_window = 128;   % +/- range around 0 lag time in which peak is sought

low_freq_cutoff = 10;
high_freq_cutoff = 100; %hz low pass cut off for integrating power

%triangle correction factor for finite trial length
triangle = (size(data.spike_data,2)*neural_bin_width/1000) - abs(time)/1000;
%spike_density_filter = [0.25 0.5 0.25];    
spike_density_filter = [1];
%trial = 0;
  
[conditions, unique_conds, param, num_conditions] = RegenerateConditionList(data, BegTrial, EndTrial, PATH, FILE, Protocol);

cross_corr_peak_time = zeros(1,num_conditions);
cross_corr_peak_power = zeros(1,num_conditions);
cross_corr_peak_height = zeros(1,num_conditions);
cross_corr_peak_area = zeros(1,num_conditions);

power_p_values = zeros(1, num_conditions);
peak_height_p_values = zeros(1, num_conditions);
peak_area_p_values = zeros(1, num_conditions);

bootstrap_peak_time = zeros(num_bootstraps, num_conditions) ;
bootstrap_peak_height = zeros(num_bootstraps, num_conditions) ;
bootstrap_peak_power = zeros(num_bootstraps, num_conditions) ;
bootstrap_peak_area = zeros(num_bootstraps, num_conditions) ;

%now, remove trials from hor_disp and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:size(conditions,2);												% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );  
  
smoothed_shuffle = zeros(num_conditions, 2*range + 1);
    
num_rows = num_conditions + 2;
num_columns = 4;
plot_offset = num_columns * 2;


output = 0;

figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [0 20 500 600], 'Name', 'Cross Correlogram Analysis');
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
line = sprintf('Window of Analysis: +/-%3d ms', peak_window);
text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Range for Averaging Peak: %2d ms', 2*peak_range + 1);
text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Freq Band For Peak Power: %2d-%2d Hz', low_freq_cutoff, high_freq_cutoff );
text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;

line = sprintf('Number of Bootstraps: %8d', num_bootstraps);
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

histo = zeros(num_conditions, size(time,2));
histo2 = zeros(num_conditions, size(time,2));
smoothed_shuffle = zeros(num_conditions, size(time,2));
 
%analyze by unique condition
for cond = 2:2 %1: num_conditions
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
      
    corr_matrix = zeros( num_reps*(num_reps), 2*range/neural_bin_width + 1);
    cross_corr_matrix = corr_matrix;
    bootstrap_corr_matrix = corr_matrix;
    
    combin = 0;
    for trial1 = 1 : num_reps
        for trial2 = 1 : num_reps
            combin = combin + 1;
            corr_matrix(combin,:) = xcorr(spikes(trial1, :), spikes2(trial2, : ), range/neural_bin_width)./triangle;		    
            %normalize by geometric mean spike rate
            if ( (neural_db1 == SPIKE_DB) | (neural_db2 == SPIKE_DB) )
               norm_factor = ( (mean( spikes(trial1, :) ) *1000  ) * (mean( spikes2(trial2, : ) )*1000 ) )^0.5;

            %normalize by minimum firing rate
%               norm_factor = min( [(mean( spikes(trial1, :) ) *1000  ) (mean( spikes2(trial2, : ) )*1000 ) ] );

                if (norm_factor ~= 0)
                    corr_matrix(combin,:) = corr_matrix(combin,:) / norm_factor;
                end    
            end    
        end
    end      

        
    % first obtain the general structure of the correlogram by creating smoothed shuffle
    % only from trial pairings that aren't cross correlograms
    crosses = [];
    crosses(linspace(1,(num_reps)^2, num_reps)) = 1;

    smoothed_shuffle(cond,:) = mean(corr_matrix(~logical(crosses) , :)  );
    histo2(cond,:) = mean(corr_matrix(logical(crosses),:) ); 

            
    %subtract smoothed shuffle
    cross_corr_matrix = corr_matrix - ones(size(corr_matrix,1),1) * smoothed_shuffle(cond,:);
    %subtract mean from single trial cross correlograms
    cross_corr_matrix = cross_corr_matrix - mean( cross_corr_matrix, 2 )*ones(1, size( cross_corr_matrix,2) );

    correlations = cross_corr_matrix(logical(crosses),:);
    for trial = 1:size(correlations,1)
        [freq2, corr_ampl] = FourierTransform_1D(0.001*[1 : (size(histo,2) - 1)], correlations(trial, 2:end), size(histo,2) - 1, 1, 0);
        corr_power{cond}(trial) = sum(corr_ampl((freq2 < high_freq_cutoff) & (freq2 > low_freq_cutoff)  ) .* corr_ampl( (freq2 < high_freq_cutoff)& (freq2 > low_freq_cutoff)   )   );
    end   
    
    
    %cross correlogram
    histo(cond,:) = mean(cross_corr_matrix(logical(crosses),:) ); 
        
    %calculate peak time and peak height
    cross_corr_peak_time(cond) = max(-1 + find(time == -peak_window) + time(find( histo(cond, find(time == -peak_window):find(time == peak_window) ) == max(histo(cond, find(time == -peak_window):find(time == peak_window) ) ) )  ) );    
    cross_corr_peak_height(cond) = histo(cond, (time == cross_corr_peak_time(cond) )  );
    
    %calculate peak power
    [freq2, ampl2] = FourierTransform_1D(0.001*[1 : (size(histo,2) - 1)], histo(cond, 2:end), size(histo,2) - 1, 1, 0);
    cross_corr_peak_power(cond) = sum(ampl2((freq2 < high_freq_cutoff) & (freq2 > low_freq_cutoff)  ) .* ampl2( (freq2 < high_freq_cutoff)& (freq2 > low_freq_cutoff)   )   );

    %calculate peak area
    peak_L = peak_range + ( find(time == cross_corr_peak_time(cond) ) <= peak_range)*( find(time == cross_corr_peak_time(cond) ) - peak_range - 1);
    peak_R = peak_range - ( find(time == cross_corr_peak_time(cond) ) > length(time) - peak_range)*( find(time == cross_corr_peak_time(cond) ) - length(time) + peak_range);
    
    peak = histo(cond, find(time == cross_corr_peak_time(cond) ) - peak_L : find(time == cross_corr_peak_time(cond) ) + peak_R  );
    cross_corr_peak_area(cond) = sum(peak(peak > 0)  );

    
    ampl5 = zeros(1, (size(corr_matrix,2) - 1)/2 );    
    bootstrap = zeros(1, size(time,2));
    
    for strap = 1:num_bootstraps
        shuffle = zeros(1, size(time,2));
        %shuffle trials   with replacement
        shuffle_index = fix(rand(1,num_reps)*combin) + 1;
        shuffle = mean(cross_corr_matrix(shuffle_index,:)  );
        
        %find peak time and peak height
        bootstrap_peak_time(strap, cond) = max( time(-1 + find(time == -peak_window) + find( shuffle(find(time == -peak_window):find(time == peak_window)) == max(shuffle(find(time == -peak_window):find(time == peak_window) ) ) ) )  );
        bootstrap_peak_height(strap,cond) = shuffle( (time == bootstrap_peak_time(strap,cond)  ) );   
       
        %bootstrap_peak_time(strap, cond)
        %calculate peak area
        peak_L = peak_range + ( find(time == bootstrap_peak_time(strap, cond) ) <= peak_range)*( find(time == bootstrap_peak_time(strap, cond) ) - peak_range - 1);
        peak_R = peak_range - ( find(time == bootstrap_peak_time(strap, cond) ) > length(time) - peak_range)*( find(time == bootstrap_peak_time(strap, cond) ) - length(time) + peak_range);

        
        peak = shuffle(find(time == bootstrap_peak_time(strap, cond) ) - peak_L : find(time == bootstrap_peak_time(strap, cond) ) + peak_R   );
        bootstrap_peak_area(strap, cond) = sum( peak(peak > 0)  );

        % calculated peak power
        [freq, ampl] = FourierTransform_1D(0.001*[1 : (length(shuffle) - 1)], shuffle(2:end), length(shuffle) - 1, 1, 0);   
        bootstrap_peak_power(strap, cond) = sum(ampl((freq < high_freq_cutoff) & (freq > low_freq_cutoff)  ) .* ampl( (freq < high_freq_cutoff)& (freq > low_freq_cutoff)   )   );

        ampl5 = ampl5 + ampl.^2;
        
        if (strap < 10)
            subplot(num_rows, num_columns, cond*num_columns - 1 + plot_offset);
            plot(freq,ampl.*ampl, 'g');
            hold on;
            subplot(num_rows, num_columns, cond*num_columns - 2 + plot_offset);
            plot(time,shuffle, 'g');
            hold on;            
        end
        
    end  

% %reassigning the identity for bootstrap shuffle correlograms    
%     for strap = 1 : num_bootstraps
%         %shuffle trials   with replacement
%         shuffle_index = fix(rand(1,combin)*combin) + 1;
%         bootstrap_corr_matrix = corr_matrix( shuffle_index, :);
%         %recalculate ensemble shuffle based on reassignment from bootstrap
%         boot_shuffle = mean(bootstrap_corr_matrix(~logical(crosses), : ) );
%         
%         %subtract smoothed shuffle from corr matrix
%         bootstrap_corr_matrix = bootstrap_corr_matrix - ones(size(bootstrap_corr_matrix,1),1) * boot_shuffle;
%         %subtract mean from single trial cross correlograms
%         bootstrap_corr_matrix = bootstrap_corr_matrix - mean( bootstrap_corr_matrix, 2 )*ones(1, size( bootstrap_corr_matrix,2) );        
%         bootstrap = mean(bootstrap_corr_matrix(logical(crosses), : ) );
%         
%         bootstrap_peak_time(strap, cond) = max( time(-1 + find(time == -peak_window) + find( bootstrap(find(time == -peak_window):find(time == peak_window)) == max(bootstrap(find(time == -peak_window):find(time == peak_window) ) ) ) )  );
%            
%         [freq, ampl] = FourierTransform_1D(0.001*[1 : (length(bootstrap) - 1)], bootstrap(2:end), length(bootstrap) - 1, 1, 0);   
%         bootstrap_peak_power(strap, cond) = sum(ampl((freq < high_freq_cutoff) & (freq > low_freq_cutoff)  ) .* ampl( (freq < high_freq_cutoff)& (freq > low_freq_cutoff)   )   );
% 
%         ampl5 = ampl5 + ampl.^2;     
%         if (strap < 10)
%             subplot(num_rows, num_columns, cond*num_columns - 1 + plot_offset);
%             plot(freq,ampl, 'g');
%             hold on;
%             subplot(num_rows, num_columns, cond*num_columns - 2 + plot_offset);
%             plot(time, bootstrap, 'g'); 
%             hold on
%         end
%         
%     end  



    ampl5 = ampl5 / num_bootstraps;
    smooth_corr_peak_power(cond) = sum(ampl5(  (freq2 < high_freq_cutoff) & (freq > low_freq_cutoff) )   );
    % do not scale cross correlogram power about 0;
 %   smooth_corr_peak_power(cond) = 0;

    cross_corr_peak_power(cond) = cross_corr_peak_power(cond) - smooth_corr_peak_power(cond);
    bootstrap_peak_power(:,cond) = bootstrap_peak_power(:,cond) - smooth_corr_peak_power(cond);    
    
    power_p_values(cond) = mean(bootstrap_peak_power(:, cond) > cross_corr_peak_power(cond) );
    peak_height_p_values(cond) = mean(bootstrap_peak_height(:,cond) > cross_corr_peak_height(cond) );
    peak_area_p_values(cond) = mean(bootstrap_peak_area(:,cond) > cross_corr_peak_area(cond) );
    
    %plot out power spectra
    subplot(num_rows, num_columns, cond*num_columns - 1 + plot_offset);
    plot(freq2,ampl2.*ampl2, 'r', freq,ampl5,'k');
    if (cond == 1)
        title('Power Spectra');
    end
    hold on;

%     subplot(num_rows, num_columns, cond*num_columns + plot_offset);
%     peaks(:,2) = bootstrap_peak_power(:,cond);
%     peaks(:,[1 3]) = NaN;
%     peaks(1:fix(num_bootstraps*0.05), 3) = cross_corr_peak_power(cond);
%     hist (peaks, 100);
%     hold on;

    subplot(num_rows, num_columns, cond*num_columns + plot_offset);
    axis([0 100 0 100]);                    
    axis('off');
    xpos = -30;
    ypos = 80;
    font_size = 8;
    bump_size = 20;
        
    line = '';
    for modality = 1: length(param) 
        line = [line param{modality}, num2str(unique_conds(cond,modality)), ' '];
    end   
    text(xpos - 60,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Peak t, Ht, A: %3d %6.3f %6.3f', cross_corr_peak_time(cond), cross_corr_peak_height(cond), cross_corr_peak_area(cond) );
    text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('P (Ht/A): %6.3f, %6.3f', peak_height_p_values(cond), peak_area_p_values(cond) );
    text(xpos,ypos, line, 'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Peak Power: %6.3f', cross_corr_peak_power(cond));
    text(xpos,ypos, line, 'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('P (Power): %6.3f', power_p_values(cond) );
    text(xpos,ypos, line, 'FontSize',font_size);		ypos = ypos - bump_size;
    
    axis('off');
    hold on;

end	% of analyzing for each condition

[1:num_conditions; cross_corr_peak_time ; cross_corr_peak_height; cross_corr_peak_area; cross_corr_peak_power; peak_height_p_values; peak_area_p_values; power_p_values]    
 

%raw_max = max( [max(histo2) max(smoothed_shuffle) ] );
%raw_min = min( [min(histo2) min(smoothed_shuffle) ] );
sub_max = 1.5 * max( max(histo) );
sub_min = 1.5 * min( min(histo) );

trial = 0;
% plotting loop
for cond = 2:2 %1: num_conditions            
    % plot smoothed shuffle and raw cross correlograms
    subplot(num_rows, num_columns, cond*num_columns - 3 + plot_offset);
    phandle = plot(time', histo2(cond,:), 'r', time', smoothed_shuffle(cond,:),  'g');
    if (cond == 1)
        title ('Raw CCG');
    end
    hold on;
    XLim([-disp_range +disp_range]);
%        	YLim([raw_min raw_max]);
    %plot Smooth Subtracted Correlograms
    subplot(num_rows, num_columns, cond*num_columns - 2 + plot_offset);
    phandle = plot(time', histo(cond,:), line_types(1));
    XLim([-disp_range +disp_range]);
    YLim([sub_min sub_max]);
    if (cond == 1)
        title('Subtracted');
    end
    hold on            
end

    %output tuning curve metrics
    if (output == 1)
        i = size(PATH,2) - 1;
        while PATH(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
            i = i - 1;
        end   
        PATHOUT = [PATH(1:i) 'Analysis\Correlograms\'];
        i = size(FILE,2) - 1;
        while FILE(i) ~='.'
            i = i - 1;
        end
        FILEOUT = [FILE(1:i) 'ccg'];
        eval(['save ' PATHOUT FILEOUT  ' cross_corr_peak_power cross_corr_peak_time cross_corr_peak_height cross_corr_peak_area bootstrap_peak_power bootstrap_peak_time bootstrap_peak_height bootstrap_peak_area power_p_values peak_height_p_values peak_area_p_values num_bootstraps range peak_range peak_window spike_density_filter low_freq_cutoff high_freq_cutoff SpikeChanA SpikeChanB start_code stop_code BegTrial EndTrial StartOffset StopOffset StartEventBin StopEventBin PATH FILE'])
    end

