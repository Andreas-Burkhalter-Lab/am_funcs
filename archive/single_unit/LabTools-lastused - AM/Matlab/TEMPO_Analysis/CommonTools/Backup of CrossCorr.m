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

peak_range = 10; % +/- range for width of averaging across peak for peak height, area, NCC
peak_window = 128;   % +/- range around 0 lag time in which peak is sought

low_freq_cutoff1 = 0;
high_freq_cutoff1 = 8; %hz low pass cut off for integrating power

low_freq_cutoff2 = 8;
high_freq_cutoff2 = 30; %hz low pass cut off for integrating power

low_freq_cutoff3 = 20;
high_freq_cutoff3 = 70; %hz low pass cut off for integrating power

%triangle correction factor for finite trial length
triangle = (size(data.spike_data,2)*neural_bin_width/1000) - abs(time)/1000;
%spike_density_filter = [1/16 4/16 6/16 4/16 1/16];    
%spike_density_filter = [0.25 0.5 0.25];    
spike_density_filter = [1];
%trial = 0;
  
[conditions, unique_conds, param, num_conditions] = RegenerateConditionList(data, BegTrial, EndTrial, PATH, FILE, Protocol);

cross_corr_peak_time = zeros(1,num_conditions);
cross_corr_peak_power = zeros(3,num_conditions);
cross_corr_peak_height = zeros(1,num_conditions);
cross_corr_peak_area = zeros(1,num_conditions);

power_p_values = zeros(1, num_conditions);
peak_height_p_values = zeros(1, num_conditions);
peak_area_p_values = zeros(1, num_conditions);

bootstrap_peak_time = zeros(num_bootstraps, num_conditions) ;
bootstrap_peak_height = zeros(num_bootstraps, num_conditions) ;
bootstrap_peak_power = zeros(num_bootstraps, num_conditions) ;
bootstrap_peak_area = zeros(num_bootstraps, num_conditions) ;
bootstrap_NCC = zeros(num_bootstraps, num_conditions) ;

%now, remove trials from hor_disp and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:size(conditions,2);												% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );  
  
smoothed_shuffle = zeros(num_conditions, 2*range + 1);
    
num_rows = num_conditions + 2;
num_columns = 4;
plot_offset = num_columns * 2;

%load msac times and regions with evoked response...
loadmsac = 1;
    if (loadmsac == 1)
        i = size(PATH,2) - 1;
        while PATH(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
            i = i - 1;
        end   
        PATHOUT = [PATH(1:i) 'Analysis\Saccades\'];
        i = size(FILE,2) - 1;
        while FILE(i) ~='.'
            i = i - 1;
        end
        FILEOUT = [FILE(1:i) 'msc'];
        eval(['load ' PATHOUT FILEOUT ' msac_end_times msac_start_times -MAT'])
        
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
        eval(['load ' PATHOUT FILEOUT ' mini_segment_start_times mini_segment_end_times -MAT'])
    
    end



output = 1;

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
line = sprintf('Freq Band For Peak Power: %2d-%2d Hz', low_freq_cutoff1, high_freq_cutoff1 );
text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Freq Band For Peak Power: %2d-%2d Hz', low_freq_cutoff2, high_freq_cutoff2 );
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
histo3 = zeros(num_conditions, size(time,2));
histo4 = zeros(num_conditions, size(time,2));
histo5 = zeros(1 , size(time,2));
histo6 = zeros(1 , size(time,2));
histo7 = zeros(1 , size(time,2));
fit_cross_correlograms = zeros(num_conditions, size(time,2));
fit_gaussian = zeros(num_conditions, size(time,2));
fit_gabor = zeros(num_conditions, size(time,2));

smoothed_shuffle = zeros(num_conditions, size(time,2));
smoothed_shuffle_GaborFit = zeros(num_conditions, size(time,2));
 





%analyze by unique condition
for cond = 2:2 %1: num_conditions
    cond
 %   for segment = 1: 1 %length(mini_segment_start_times{2})
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
    
if ( (neural_db1 == LFP_DB) | (neural_db2 == LFP_DB) )
    c_coef = zeros(num_reps, 2 * range/neural_bin_width  + 1);
    for trial1 = 1:num_reps
        for tau = -range/neural_bin_width : 0
            shifted1 = spikes(trial1, 1 : end + tau);
            shifted2 = spikes2(trial1, 1 - tau : end);         
            cc_matrix = corrcoef(shifted1, shifted2  );
            c_coef(trial1, tau + range/neural_bin_width + 1) = cc_matrix(1,2);
        end    

        for tau = 1 : range/neural_bin_width
            shifted1 = spikes(trial1, 1 + tau: end);
            shifted2 = spikes2(trial1, 1 : end - tau);
            cc_matrix = corrcoef(shifted1, shifted2  );
            c_coef(trial1, tau + range/neural_bin_width + 1) = cc_matrix(1,2);
        end    

    end
    
    histo(cond,:) = mean(c_coef ); 

        %calculate peak power
    [freq2, ampl2] = FourierTransform_1D( (neural_bin_width/1000 )*[1 : (size(histo,2) - 1)], histo(cond, 2:end), size(histo,2) - 1, 1, 0);

    cross_corr_peak_power(1, cond) = sum(ampl2((freq2 < high_freq_cutoff1) & (freq2 > low_freq_cutoff1)  ) .* ampl2( (freq2 < high_freq_cutoff1)& (freq2 > low_freq_cutoff1)   )   );

    cross_corr_peak_power(2, cond) = sum(ampl2((freq2 < high_freq_cutoff2) & (freq2 > low_freq_cutoff2)  ) .* ampl2( (freq2 < high_freq_cutoff2)& (freq2 > low_freq_cutoff2)   )   );

    cross_corr_peak_power(3, cond) = sum(ampl2((freq2 < high_freq_cutoff3) & (freq2 > low_freq_cutoff3)  ) .* ampl2( (freq2 < high_freq_cutoff3)& (freq2 > low_freq_cutoff3)   )   );
    
    cross_corr_peak_time(cond) = max(-1 + find(time == -peak_window) + time(find( histo(cond, find(time == -peak_window):find(time == peak_window) ) == max(histo(cond, find(time == -peak_window):find(time == peak_window) ) ) )  ) );    
    cross_corr_peak_height(cond) = histo(cond, (time == cross_corr_peak_time(cond) )  );

    
else    %analyze spike DB
    corr_matrix = zeros( num_reps*(num_reps), 2*range/neural_bin_width + 1);
    auto_corr_matrix1 = zeros (num_reps*num_reps, 2*range/neural_bin_width + 1);
    auto_corr_matrix2 = zeros (num_reps*num_reps, 2*range/neural_bin_width + 1);
    
    cross_corr_matrix = corr_matrix;
    corr_matrix_GaborFit = corr_matrix;
    bootstrap_corr_matrix = corr_matrix;
    
    combin = 0;
    for trial1 = 1 : num_reps
        for trial2 = 1 : num_reps
            combin = combin + 1;
            corr_matrix(combin,:) = xcorr(spikes(trial1, :), spikes2(trial2, : ), range/neural_bin_width); 
            %convert to coincidences / sec by normalizing for finite trial length
            corr_matrix_GaborFit(combin,:) = corr_matrix(combin,:);
            corr_matrix(combin,:) = corr_matrix(combin,:)./ triangle;		    
            auto_corr_matrix1(combin, :) = (  xcorr(spikes(trial1, :), spikes(trial2, : ), range/neural_bin_width)   ) ./ triangle  ;		    
            auto_corr_matrix2(combin, :) = (  xcorr(spikes2(trial1, :), spikes2(trial2, : ), range/neural_bin_width)  )./ triangle  ;		    
          
            if ( (neural_db1 == SPIKE_DB) | (neural_db2 == SPIKE_DB) )
                %normalize by geometric mean spike rate
                norm_factor = ( (mean( spikes(trial1, :) ) *1000  ) * (mean( spikes2(trial2, : ) )*1000 ) )^0.5;
                %generate autocorrelogram matrix
                auto_corr_matrix1(combin, :) = ( auto_corr_matrix1(combin, :) ) / (mean( spikes(trial1, :) ) *1000  );  
                auto_corr_matrix2(combin, :) = ( auto_corr_matrix2(combin, :)  ) / (mean( spikes2(trial2, : ) )*1000 );
                %normalize by minimum firing rate
%                norm_factor = min( [(mean( spikes(trial1, :) ) *1000  ) (mean( spikes2(trial2, : ) )*1000 ) ] );
                
 %               norm_factor = 1;
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

    %histo 2 = raw cross correlogram no shuffle subtraction
    smoothed_shuffle(cond,:) = mean(corr_matrix(~logical(crosses) , :)  );
    histo2(cond,:) = mean(corr_matrix(logical(crosses),:) ); 

    % calculate mean autocorrelogram and subtract shuffle1
    smoothed_auto_shuffle1(cond,:) = mean(auto_corr_matrix1(~logical(crosses), :) );
    smoothed_auto_shuffle2(cond,:) = mean(auto_corr_matrix2(~logical(crosses), :) );
    
   auto_corr_matrix1 = auto_corr_matrix1 - ones(size( auto_corr_matrix1,1), 1  )*smoothed_auto_shuffle1(cond,:);
   auto_corr_matrix2 = auto_corr_matrix2 - ones(size( auto_corr_matrix2,1), 1  )*smoothed_auto_shuffle2(cond,:);
           
    histo3(cond,:) = mean(auto_corr_matrix1(logical(crosses), :) );
    histo4(cond,:) = mean(auto_corr_matrix2(logical(crosses), :) );
    
    auto_corr_area(cond,1) = sum( histo3(cond, find(time == 0 ) - peak_range / neural_bin_width : find(time == 0 ) + peak_range / neural_bin_width )  );
    auto_corr_area(cond,2) = sum( histo4(cond, find(time == 0 ) - peak_range / neural_bin_width : find(time == 0 ) + peak_range / neural_bin_width )  );

    %subtract smoothed shuffle
    cross_corr_matrix = corr_matrix - ones(size(corr_matrix,1),1) * smoothed_shuffle(cond,:);
    %subtract mean from single trial cross correlograms
    cross_corr_matrix = cross_corr_matrix - mean( cross_corr_matrix, 2 )*ones(1, size( cross_corr_matrix,2) );

  
    %cross correlogram = mean of the subtracted cross correlograms
    histo(cond,:) = mean(cross_corr_matrix(logical(crosses),:) ); 

    smoothed_shuffle_GaborFit(cond,:) = sum(corr_matrix_GaborFit(~logical(crosses) , :)  ) / num_reps;
    % histo7 is raw cross correlograms no triangle but with shuffle
    % subtraction (needed because shuffle not flat) but with mean of
    % shuffle added
    histo7(cond,:) = sum(corr_matrix_GaborFit(logical(crosses),:) ) - smoothed_shuffle_GaborFit(cond,:) + mean(smoothed_shuffle_GaborFit(cond,:) ); 
        
    %calculate peak time and peak height
    cross_corr_peak_time(cond) = max(-1 + find(time == -peak_window) + time(find( histo(cond, find(time == -peak_window):find(time == peak_window) ) == max(histo(cond, find(time == -peak_window):find(time == peak_window) ) ) )  ) );    
    cross_corr_peak_height(cond) = histo(cond, (time == cross_corr_peak_time(cond) )  );
    
    %calculate peak power
    [freq2, ampl2] = FourierTransform_1D( (neural_bin_width/1000 )*[1 : (size(histo,2) - 1)], histo(cond, 2:end), size(histo,2) - 1, 1, 0);

    if cond == 2
        power_single_object = ampl2.*ampl2;
    end
    
    cross_corr_peak_power(1, cond) = sum(ampl2((freq2 < high_freq_cutoff1) & (freq2 > low_freq_cutoff1)  ) .* ampl2( (freq2 < high_freq_cutoff1)& (freq2 > low_freq_cutoff1)   )   );

    cross_corr_peak_power(2, cond) = sum(ampl2((freq2 < high_freq_cutoff2) & (freq2 > low_freq_cutoff2)  ) .* ampl2( (freq2 < high_freq_cutoff2)& (freq2 > low_freq_cutoff2)   )   );

    cross_corr_peak_power(3, cond) = sum(ampl2((freq2 < high_freq_cutoff3) & (freq2 > low_freq_cutoff3)  ) .* ampl2( (freq2 < high_freq_cutoff3)& (freq2 > low_freq_cutoff3)   )   );
    
    %calculate peak area around peak height
%    peak_L = peak_range + ( find(time == cross_corr_peak_time(cond) ) <= peak_range)*( find(time == cross_corr_peak_time(cond) ) - peak_range - 1);
%    peak_R = peak_range - ( find(time == cross_corr_peak_time(cond) ) > length(time) - peak_range)*( find(time == cross_corr_peak_time(cond) ) - length(time) + peak_range);    
%    peak = histo(cond, find(time == cross_corr_peak_time(cond) ) - peak_L : find(time == cross_corr_peak_time(cond) ) + peak_R  );

    %calculate peak area around 0 ms time lag
    cross_corr_peak_area(cond) = sum(  histo(cond, find(time == 0 ) - peak_range/neural_bin_width : find(time == 0 ) + peak_range/neural_bin_width  )    );

    NCC(cond) = cross_corr_peak_area(cond)/(auto_corr_area(cond,1) * auto_corr_area(cond,2) )^0.5;
    
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

        histo5 = mean(auto_corr_matrix1( shuffle_index, :) );
        histo6 = mean(auto_corr_matrix2( shuffle_index, :) );
           
        bootstrap_NCC(strap, cond) = bootstrap_peak_area(strap, cond)/( sum( histo6(find(time == 0 ) - peak_range : find(time == 0 ) + peak_range )  ) *         sum( histo5(find(time == 0 ) - peak_range : find(time == 0 ) + peak_range )  )   )^0.5;
        
        % calculated peak power
        [freq, ampl] = FourierTransform_1D( (neural_bin_width/1000 ) *[1 : (length(shuffle) - 1)], shuffle(2:end), length(shuffle) - 1, 1, 0);   
        bootstrap_peak_power(strap, cond) = sum(ampl((freq < high_freq_cutoff1) & (freq > low_freq_cutoff1)  ) .* ampl( (freq < high_freq_cutoff1)& (freq > low_freq_cutoff1)   )   );

        ampl5 = ampl5 + ampl.^2;
        
        if (strap < 10)
 %           subplot(num_rows, num_columns, cond*num_columns - 1 + plot_offset);
 %           plot(freq,ampl.*ampl, 'g');
 %           hold on;
            subplot(num_rows, num_columns, cond*num_columns - 3 + plot_offset);
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
%     smooth_corr_peak_power(cond) = sum(ampl5(  (freq2 < high_freq_cutoff1) & (freq2 > low_freq_cutoff1) )   );
%     % do not scale cross correlogram power about 0;
% %    smooth_corr_peak_power(cond) = 0;
% 
%     cross_corr_peak_power(cond) = cross_corr_peak_power(cond) - smooth_corr_peak_power(cond);
%     bootstrap_peak_power(:,cond) = bootstrap_peak_power(:,cond) - smooth_corr_peak_power(cond);    
%     

    power_p_values(cond) = mean(bootstrap_peak_power(:, cond) > cross_corr_peak_power(cond) );
    peak_height_p_values(cond) = mean(bootstrap_peak_height(:,cond) > cross_corr_peak_height(cond) );
    peak_area_p_values(cond) = mean(bootstrap_peak_area(:,cond) > cross_corr_peak_area(cond) );
    
    NCC_p_values(cond) = mean(bootstrap_NCC(:,cond) > NCC(cond) );
    
    %plot out power spectra
%    subplot(num_rows, num_columns, cond*num_columns - 1 + plot_offset);
%    plot(freq2,ampl2.*ampl2, 'r', freq2,ampl5,'k');
%    xlim( [0 max(freq2)] );
%    xlim( [0 100] );
%    if (cond == 1)
%        title('Power Spectra');
%    end
%    hold on;

    gabor_fit = 0;
    NC(cond) = NaN;
    pars{cond} = NaN;
    raw = [];
    fit_cross_correlogram(cond,:) = 0*time;
    if (gabor_fit == 1)
        raw(:,1) = reshape([time' *ones(1 , num_reps)], 1, num_reps * length(time) )' ;
        raw(:,2) = reshape(corr_matrix_GaborFit(logical(crosses),:)', 1, num_reps * length(time))' ;

        fixed_param_flags = zeros(8,1); %by default, all 8 parameters will vary
        fixed_param_values = zeros(8,1); %override these values and flags to fix a parameter
        [pars{cond},freq(cond)] = gaborgauss_sumfit([time' histo7(cond,:)'],raw,fixed_param_flags,fixed_param_values);
        fit_cross_correlogram(cond,:) = real (gaborgauss_sumfunc(time, pars{cond}) ) ;
        if pars{cond}(6) <= 1/pars{cond}(5)
            A(cond) = pars{cond}(2) + pars{cond}(7);
            O(cond) = pars{cond}(1);
        else
            A(cond) = pars{cond}(2);
            O(cond) = pars{cond}(1) + pars{cond}(7);
        end
        NC(cond) = A(cond)/O(cond) * 100;
        fit_gabor(cond,:) = pars{cond}(1) + pars{cond}(2)*exp( -(  abs(time - pars{cond}(3)) / pars{cond}(4)).^pars{cond}(8)).*cos(2*pi*pars{cond}(5)/1000 * (time - pars{cond}(3)) + 2*pi/1000 * pars{cond}(6) ); 
        fit_gaussian(cond,:) = pars{cond}(1) + pars{cond}(7) * exp( -((time - pars{cond}(3) )/ pars{cond}(9) ).^2);

    end
    
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
    line = sprintf('Peak Power: %6.3f %6.3f %6.3f', cross_corr_peak_power(1, cond), cross_corr_peak_power(2, cond), cross_corr_peak_power(3, cond));
    text(xpos,ypos, line, 'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('P (Power): %6.3f', power_p_values(cond) );
    text(xpos,ypos, line, 'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('NCC: %6.3f   NC: %6.3f', NCC(cond), NC(cond) );
    text(xpos,ypos, line, 'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Curve Fit Params: %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f', pars{cond} );
    text(xpos - 60,ypos, line, 'FontSize',font_size);		ypos = ypos - bump_size;
    
    axis('off');
    hold on;
    
    
 
end %SpikeDB

end	% of analyzing for each condition

% [1:num_conditions; cross_corr_peak_time ; cross_corr_peak_height; cross_corr_peak_area; cross_corr_peak_power(1,:); cross_corr_peak_power(2,:); peak_height_p_values; peak_area_p_values; power_p_values]    


if ( (neural_db1 == LFP_DB) | (neural_db2 == LFP_DB) )

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
        eval(['save ' PATHOUT FILEOUT  ' freq2 cross_corr_peak_power cross_corr_peak_time cross_corr_peak_height range peak_range peak_window spike_density_filter low_freq_cutoff1 high_freq_cutoff1 low_freq_cutoff2 high_freq_cutoff2 SpikeChanA SpikeChanB start_code stop_code BegTrial EndTrial StartOffset StopOffset StartEventBin StopEventBin PATH FILE'])
    end

    
else    %spike DB    

%raw_max = max( [max(histo2) max(smoothed_shuffle) ] );
%raw_min = min( [min(histo2) min(smoothed_shuffle) ] );
sub_max = 1.5 * max( max(histo) );
sub_min = 1.5 * min( min(histo) );

trial = 0;
% plotting loop
for cond = 2:2 %1: num_conditions            
%    % plot smoothed shuffle and raw cross correlograms
%    subplot(num_rows, num_columns, cond*num_columns - 3 + plot_offset);
%    phandle = plot(time', smoothed_shuffle(cond,:),  'g', time', histo2(cond,:), 'r');
%    if (cond == 1)
%        title ('Raw CCG');
%    end
%    hold on;
    XLim([-disp_range +disp_range]);
    subplot(num_rows, num_columns, cond*num_columns - 3 + plot_offset);
    phandle = plot(time', histo(cond,:), 'r');

    XLim([-disp_range +disp_range]);
    YLim([sub_min sub_max]);
    if (cond == 1)
        title('Subtracted');
    end
    hold on         
    
    % plot smoothed shuffle and raw cross correlograms
    subplot(num_rows, num_columns, cond*num_columns - 2 + plot_offset);
    phandle = plot(time', histo7(cond,:), 'r', time, fit_cross_correlogram(cond,:), 'g');
    XLim([-disp_range +disp_range]);
 %   YLim([sub_min sub_max]);

    if (cond == 1)
        title ('Unnormalized');
    end
    hold on; 
    
   % plot gabor and gaussian parts of fit correlogram
    subplot(num_rows, num_columns, cond*num_columns - 1 + plot_offset);
    
    phandle = plot(time', fit_gabor(cond,:) , 'r', time, fit_gaussian(cond,:) , 'g');
    XLim([-disp_range +disp_range]);

    if (cond == 1)
        title('Fit Correlogram Components');
    end

    hold on;
    
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
        eval(['save ' PATHOUT FILEOUT  ' freq2 power_single_object cross_corr_peak_power cross_corr_peak_time cross_corr_peak_height cross_corr_peak_area NCC NC bootstrap_peak_power bootstrap_peak_time bootstrap_peak_height bootstrap_peak_area power_p_values peak_height_p_values peak_area_p_values num_bootstraps range peak_range peak_window spike_density_filter low_freq_cutoff1 high_freq_cutoff1 low_freq_cutoff2 high_freq_cutoff2 SpikeChanA SpikeChanB start_code stop_code BegTrial EndTrial StartOffset StopOffset StartEventBin StopEventBin PATH FILE'])
    end

end