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

%regenerate condition list
[conditions, unique_conds, param, num_conditions] = RegenerateConditionList(data, BegTrial, EndTrial, PATH, FILE, Protocol);

eval (['warning off MATLAB:divideByZero']);

if Protocol > 200
    %load msac times and regions with evoked response...
    loadmsac =1;
    switch Protocol
        case FIX_1_23_45
            comp_corr_conds = [2 5 8];
            num_comp=2;
        case FIX_1_23
            comp_corr_conds = [2 5];
            num_comp=1;
        case FIX_VARY_HISTORY
            comp_corr_conds = [3 4 7 8];
            num_comp=2;
        case FIX_VARY_BACKCOLOR    
            comp_corr_conds = [2 3];
            num_comp=1;
    end  
else
    loadmsac = 0;
end

if (loadmsac == 1)
    %load microsaccades
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
    
    %load PSTH
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
    eval(['load ' PATHOUT FILEOUT ' selected_spikes1 selected_spikes2 selected_FR1 selected_FR2 time hist_data2 hist_data2B minimum_segment_length mini_segment_start_times mini_segment_end_times -MAT'])
    plot_PSTH = 1;
    hist_time = time;
    time = [];
else
    selected_spikes1 = [];
    selected_spikes2 = [];
    selected_FR1 = [];
    selected_FR2 = [];
    minimum_segment_length = min(StopEventBin - StartEventBin);
    for cond = 1:num_conditions
        mini_segment_start_times{cond} = 0;
        mini_segment_end_times{cond} = min(StopEventBin - StartEventBin);
    end        
    plot_PSTH = 0;
end

%keep track of spike channels before reindexing
SpikeChanA = SpikeChan;
SpikeChanB = SpikeChan2;
num_spike_channels = size(data.spike_data,1);

%needs to be revised
if Protocol > 200 %change spike DB
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
else
    neural_db1 = SPIKE_DB;
    neural_db2 = SPIKE_DB;   
end



%first, we'll need the bin_width (in sec) of our spike raster + events log
h = data.htb_header{neural_db1};	%for convenience
neural_bin_width = (h.skip + 1) / (h.speed_units / h.speed) * 1000;  %bin width of neural signal in ms

AnalInterval = min( (StopEventBin + StopOffset) - (StartEventBin+ StartOffset) )  ;

rand('state',sum(100*clock));

num_permutations = 1000;
num_permutations2 = 1000;

range = 50;		% Range for generating cross correlograms
disp_range = 50;    %Range for histograms
time = -range:neural_bin_width:+range; %actual time range for cross correlograms
correlogram_norm_type = 1; %0 = no norm, 1 = GMSR, 2 = min firing rate of the two cells, 3 product of firing rates;

peak_range = 5; % +/- range for width of averaging across peak for peak height, area, NCC
peak_range2 = 10; % +/- range for width of averaging across peak for peak height, area, NCC
peak_range3 = 20; % +/- range for width of averaging across peak for peak height, area, NCC
peak_range4 = 40; % +/- range for width of averaging across peak for peak height, area, NCC
peak_range5 = 2; % +/- range for width of averaging across peak for peak height, area, NCC

peak_window = 10;   % +/- range around 0 lag time in which peak is sought

%spike_density_filter = [1/16 4/16 6/16 4/16 1/16];    
%spike_density_filter = [0.25 0.5 0.25];    
spike_density_filter = [1];

cross_corr_peak_time = zeros(1,num_conditions);
cross_corr_peak_height = zeros(1,num_conditions);
cross_corr_peak_area = zeros(1,num_conditions);

NCC = zeros(1, num_conditions);
NCC2 = zeros(1, num_conditions);
NCC3 = zeros(1, num_conditions);
NCC4 = zeros(1, num_conditions);
NCC5 = zeros(1, num_conditions);

NC = zeros(1, num_conditions);

permutation_CCG = zeros(num_permutations, size(time,2));
compare_corr_CCG1 = zeros(num_permutations2, size(time,2));
compare_corr_CCG2 = zeros(num_permutations2, size(time,2));


permutation_CCG_area = zeros(num_permutations,1);
permutation_ACG1_area = zeros(num_permutations,1);
permutation_ACG2_area = zeros(num_permutations,1);
permutation_NCC = zeros(num_permutations,1);
permutation_NCC_p_values = zeros(1, num_conditions);
permutation_NCC2_p_values = zeros(1, num_conditions);
permutation_NCC3_p_values = zeros(1, num_conditions);
permutation_NCC4_p_values = zeros(1, num_conditions);
permutation_NCC5_p_values = zeros(1, num_conditions);

Raw_CCG = zeros(num_conditions, size(time,2));
Accounted_Var = zeros(1,num_conditions);

%now, remove trials from hor_disp and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:size(conditions,2);												% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );  

smoothed_shuffle = zeros(num_conditions, 2*range + 1);

num_rows = num_conditions + 2;
num_columns = 4;
plot_offset = num_columns * 2;

num_selected_segments = zeros(1,num_conditions);
num_reps = zeros(1,num_conditions);

output = 1;
gabor_fit = 1;
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
    if exist([PATHOUT FILEOUT])
        eval(['load ' PATHOUT FILEOUT ' -MAT'])
    end    
end  
curr_fig = figure;
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
line = sprintf('Range for Averaging Peak: %2d ms', 2*peak_range2 + 1);
text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Range for Averaging Peak: %2d ms', 2*peak_range3 + 1);
text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Range for Averaging Peak: %2d ms', 2*peak_range4 + 1);
text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Range for Averaging Peak: %2d ms', 2*peak_range5 + 1);
text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;

line = sprintf('Number of permutations: %8d', num_permutations);
text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;

xpos = 40;
ypos = 150;
line = sprintf('Protocol: %d (%s)', Protocol, protocol_names{Protocol+1});
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Trial range for analysis: %d -> %d', BegTrial, EndTrial);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Spike channels analyzed: %d, %d', SpikeChan, SpikeChan2);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Analysis Start Event: %s, offset %d ms', event_names{start_code}, StartOffset);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;
line = sprintf('Analysis Stop Event: %s, offset %d ms', event_names{stop_code}, StopOffset);
text(xpos,ypos,line,'FontSize',font_size);		ypos = ypos - bump_size;

histo = zeros(num_conditions, size(time,2));
histo3 = zeros(num_conditions, size(time,2));
histo4 = zeros(num_conditions, size(time,2));
histo5 = zeros(1 , size(time,2));
histo6 = zeros(1 , size(time,2));
histo7 = zeros(1 , size(time,2));
ACG1 = histo;
ACG2 = histo;
CCG = histo;
fit_cross_correlograms = zeros(num_conditions, size(time,2));
fit_gaussian = zeros(num_conditions, size(time,2));
fit_gabor = zeros(num_conditions, size(time,2));

smoothed_shuffle = zeros(num_conditions, size(time,2));
smoothed_auto_shuffle1 = zeros(num_conditions, size(time,2));
smoothed_auto_shuffle2 = zeros(num_conditions, size(time,2));

smoothed_shuffle_GaborFit = zeros(num_conditions, size(time,2));

% segment_start_times = 1:step_increment:3000;
% segment_end_times = segment_start_times + segment_length - 1; 

%For comparing corr strength, need to determine number of reps in each
%conditions
for cond =1 :num_conditions
    cond
    pars{cond} = NaN;
    
    
    SetTrials(1,:) = (conditions(1,:) == unique_conds(cond,1) & select_trials); 	    
    for modality = 1: size(conditions,1) - 1         	
        NextSetTrials = (conditions(modality,:) == unique_conds(cond,modality) & select_trials); 
        SetTrials = SetTrials & NextSetTrials;				
    end	 
    
    repetitions{cond} = find(SetTrials==1); 
    num_reps(cond) = length(repetitions{cond});	
end

shuffle_index2 = zeros(num_permutations2,num_reps(comp_corr_conds(1))  + num_reps(comp_corr_conds(2))           );
shuffle_index3 = zeros(num_permutations2,num_reps(comp_corr_conds(1)) + num_reps(comp_corr_conds(2))              );

for cond = 1:length (comp_corr_conds)
    comp_corr_CCG{cond } = zeros(num_reps(comp_corr_conds(cond )), size(time,2) );
    comp_corr_ACG1{cond } = zeros(num_reps(comp_corr_conds(cond )), size(time,2) );
    comp_corr_ACG2{cond } = zeros(num_reps(comp_corr_conds(cond )), size(time,2) );
    
end

%analyze by unique condition
for cond = 1 :num_conditions
    cond
    pars{cond} = NaN;
    
   
    reps = repetitions{cond};
    
    num_selected_segments(cond) = length(mini_segment_start_times{cond});
    CCG_segment = zeros(num_selected_segments(cond), length(time)); 
    if (num_selected_segments(cond) == 1)
        mini_segment_start_times{cond} = 1:100:mini_segment_end_times{cond}-100;
        mini_segment_end_times{cond} = mini_segment_start_times{cond} + 100;
    end
    
    shuffle_index{cond} = zeros(num_permutations,num_reps(cond));
    %first set up permutations
    for perm = 1 : num_permutations
        %shuffle trials   with replacement
        shuffle_index{cond}(perm,:) = 1:num_reps(cond):num_reps(cond)*num_reps(cond);
        shuffle_index{cond}(perm,:) = randperm(num_reps(cond)) - 1 + shuffle_index{cond}(perm,:);
    end  

   
   for segment = 1: length(mini_segment_start_times{cond} )
        % These next few lines check for trials belonging to a single condition
        segment
        
        spikes = [];    
        spikes2 = [];
        % Now do cross correlograms 
        for trial = 1:num_reps(cond)
            start_spikebin = floor( StartEventBin(reps(trial)) / neural_bin_width ) + floor ( mini_segment_start_times{cond}(segment)  / neural_bin_width );
            stop_spikebin = floor( StartEventBin(reps(trial)) / neural_bin_width ) + floor ( mini_segment_end_times{cond}(segment) / neural_bin_width  );
            
            switch neural_db1
                case SPIKE_DB
                    spikes(trial,:) = conv(spike_density_filter, data.spike_data(SpikeChan, start_spikebin: stop_spikebin, reps(trial))  ); 
                case LFP_DB
                    spikes(trial,:) = data.lfp_data(SpikeChan, start_spikebin: stop_spikebin  , reps(trial)); 
            end   
            switch neural_db2
                case SPIKE_DB
                    spikes2(trial,:) = conv(spike_density_filter, data.spike_data(SpikeChan2, start_spikebin: stop_spikebin, reps(trial))  );     
                case LFP_DB
                    spikes2(trial,:) = data.lfp_data(SpikeChan2, start_spikebin: stop_spikebin   , reps(trial)); 
            end   
        end   
        
        corr_matrix = zeros( num_reps(cond)*(num_reps(cond)), 2*range/neural_bin_width + 1);
        auto_corr_matrix1 = zeros (num_reps(cond)*num_reps(cond), 2*range/neural_bin_width + 1);
        auto_corr_matrix2 = zeros (num_reps(cond)*num_reps(cond), 2*range/neural_bin_width + 1);
               
        [histo(cond,:) smoothed_shuffle(cond,:) corr_matrix CCG_segment(segment,:)] = compute_correlogram(spikes, spikes2, time, neural_bin_width, correlogram_norm_type);
        [histo3(cond,:) smoothed_auto_shuffle1(cond,:) auto_corr_matrix1] = compute_correlogram(spikes, spikes, time, neural_bin_width, correlogram_norm_type);
        [histo4(cond,:) smoothed_auto_shuffle2(cond,:) auto_corr_matrix2] = compute_correlogram(spikes2, spikes2, time, neural_bin_width, correlogram_norm_type);
                
        %calculate peak area around 0 ms time lag
        cross_corr_peak_area(cond) = sum(  histo(cond, find(time == 0 ) - peak_range/neural_bin_width : find(time == 0 ) + peak_range/neural_bin_width  )    );
        
        %For detecting significance of synchrony in a single data set
        for perm = 1:num_permutations
            permutation_CCG(perm,:) = permutation_CCG(perm,:) + mean(corr_matrix(shuffle_index{cond}(perm,:),:)  );
        end  
        
        if ~isempty(find(cond == comp_corr_conds) )
                comp_corr_CCG{find(cond == comp_corr_conds)} = comp_corr_CCG{find(cond == comp_corr_conds)} + corr_matrix(linspace(1,(num_reps(cond))^2, num_reps(cond)), :);
                comp_corr_ACG1{find(cond == comp_corr_conds)} = comp_corr_ACG1{find(cond == comp_corr_conds)} + auto_corr_matrix1(linspace(1,(num_reps(cond))^2, num_reps(cond)), :);
                comp_corr_ACG2{find(cond == comp_corr_conds)} = comp_corr_ACG2{find(cond == comp_corr_conds)} + auto_corr_matrix2(linspace(1,(num_reps(cond))^2, num_reps(cond)), :);
        end
        
        ACG1(cond, :) = ACG1(cond, :) + histo3(cond,:);
        ACG2(cond, :) = ACG2(cond, :) + histo4(cond,:);
        CCG(cond, :) =  CCG(cond, :) + histo(cond,:);
        
    end %end mini-segment
    
    Raw_CCG(cond,:) = sum (CCG_segment, 1);
    ACG1(cond, :) = ACG1(cond, :) / segment;
    ACG2(cond, :) = ACG2(cond, :) / segment;
    CCG(cond, :) =  CCG(cond, :) / segment;
    permutation_CCG = permutation_CCG/segment; 
    
    auto_corr_area(cond) = sum( ACG1(cond , find(time == 0 ) - peak_range5 / neural_bin_width : find(time == 0 ) + peak_range5 / neural_bin_width )  );
    auto_corr_area2(cond) = sum( ACG2(cond , find(time == 0 ) - peak_range5 / neural_bin_width : find(time == 0 ) + peak_range5 / neural_bin_width )  );
    cross_corr_area(cond) = sum( CCG(cond , find(time == 0 ) - peak_range5 / neural_bin_width : find(time == 0 ) + peak_range5 / neural_bin_width )  );
    NCC5(cond) = cross_corr_area(cond) / (auto_corr_area(cond) * auto_corr_area2(cond) ).^0.5;

    permutation_CCG_area = sum( permutation_CCG(: , find(time == 0 ) - peak_range5 / neural_bin_width : find(time == 0 ) + peak_range5 / neural_bin_width ), 2  );
    permutation_NCC = permutation_CCG_area ./ (auto_corr_area(cond).*auto_corr_area2(cond)).^0.5;
    permutation_NCC5_p_values(cond) = sum(NCC5(cond) < permutation_NCC )/num_permutations;
    
    auto_corr_area(cond) = sum( ACG1(cond , find(time == 0 ) - peak_range4 / neural_bin_width : find(time == 0 ) + peak_range4 / neural_bin_width )  );
    auto_corr_area2(cond) = sum( ACG2(cond , find(time == 0 ) - peak_range4 / neural_bin_width : find(time == 0 ) + peak_range4 / neural_bin_width )  );
    cross_corr_area(cond) = sum( CCG(cond , find(time == 0 ) - peak_range4 / neural_bin_width : find(time == 0 ) + peak_range4 / neural_bin_width )  );
    NCC4(cond) = cross_corr_area(cond) / (auto_corr_area(cond) * auto_corr_area2(cond) ).^0.5;

    permutation_CCG_area = sum( permutation_CCG(: , find(time == 0 ) - peak_range4 / neural_bin_width : find(time == 0 ) + peak_range4 / neural_bin_width ), 2  );
    permutation_NCC = permutation_CCG_area ./ (auto_corr_area(cond).*auto_corr_area2(cond)).^0.5;
    permutation_NCC4_p_values(cond) = sum(NCC4(cond) < permutation_NCC )/num_permutations;
 
    auto_corr_area(cond) = sum( ACG1(cond , find(time == 0 ) - peak_range3 / neural_bin_width : find(time == 0 ) + peak_range3 / neural_bin_width )  );
    auto_corr_area2(cond) = sum( ACG2(cond , find(time == 0 ) - peak_range3 / neural_bin_width : find(time == 0 ) + peak_range3 / neural_bin_width )  );
    cross_corr_area(cond) = sum( CCG(cond , find(time == 0 ) - peak_range3 / neural_bin_width : find(time == 0 ) + peak_range3 / neural_bin_width )  );        
    NCC3(cond) = cross_corr_area(cond) / (auto_corr_area(cond) * auto_corr_area2(cond) ).^0.5;

    permutation_CCG_area = sum( permutation_CCG(: , find(time == 0 ) - peak_range3 / neural_bin_width : find(time == 0 ) + peak_range3 / neural_bin_width ), 2  );
    permutation_NCC = permutation_CCG_area ./ (auto_corr_area(cond).*auto_corr_area2(cond)).^0.5;
    permutation_NCC3_p_values(cond) = sum(NCC3(cond) < permutation_NCC )/num_permutations;
   
    
    auto_corr_area(cond) = sum( ACG1(cond , find(time == 0 ) - peak_range2 / neural_bin_width : find(time == 0 ) + peak_range2 / neural_bin_width )  );
    auto_corr_area2(cond) = sum( ACG2(cond , find(time == 0 ) - peak_range2 / neural_bin_width : find(time == 0 ) + peak_range2 / neural_bin_width )  );
    cross_corr_area(cond) = sum( CCG(cond , find(time == 0 ) - peak_range2 / neural_bin_width : find(time == 0 ) + peak_range2 / neural_bin_width )  );
    NCC2(cond) = cross_corr_area(cond) / (auto_corr_area(cond) * auto_corr_area2(cond) ).^0.5;

    permutation_CCG_area = sum( permutation_CCG(: , find(time == 0 ) - peak_range2 / neural_bin_width : find(time == 0 ) + peak_range2 / neural_bin_width ), 2  );
    permutation_NCC = permutation_CCG_area ./ (auto_corr_area(cond).*auto_corr_area2(cond)).^0.5;
    permutation_NCC2_p_values(cond) = sum(NCC2(cond) < permutation_NCC )/num_permutations;

    auto_corr_area(cond) = sum( ACG1(cond , find(time == 0 ) - peak_range / neural_bin_width : find(time == 0 ) + peak_range / neural_bin_width )  );
    auto_corr_area2(cond) = sum( ACG2(cond , find(time == 0 ) - peak_range / neural_bin_width : find(time == 0 ) + peak_range / neural_bin_width )  );
    cross_corr_area(cond) = sum( CCG(cond , find(time == 0 ) - peak_range / neural_bin_width : find(time == 0 ) + peak_range / neural_bin_width )  );
    NCC(cond) = cross_corr_area(cond) / (auto_corr_area(cond) * auto_corr_area2(cond) ).^0.5;
      
    %calculate peak time and peak height
    cross_corr_peak_time(cond) = max(-1 + find(time == -peak_window) + time(find( CCG(cond, find(time == -peak_window):find(time == peak_window) ) == max(CCG(cond, find(time == -peak_window):find(time == peak_window) ) ) )  ) );    
    cross_corr_peak_height(cond) = CCG(cond, (time == cross_corr_peak_time(cond) )  );
    
    permutation_CCG_area = sum( permutation_CCG(: , find(time == 0 ) - peak_range / neural_bin_width : find(time == 0 ) + peak_range / neural_bin_width ), 2  );
    permutation_NCC = permutation_CCG_area ./ (auto_corr_area(cond).*auto_corr_area2(cond)).^0.5;
    permutation_NCC_p_values(cond) = sum(NCC(cond) < permutation_NCC )/num_permutations;
    
    %now that we have cross correlogram, fit with gabor and gaussian a la
    %Konig, 1994
    NC(cond) = NaN;
    pars{cond} = NaN;
    raw = [];
    fit_cross_correlogram(cond,:) = 0*time;
    if (gabor_fit == 1)
        raw = 0;
        %fit-correlogram differences weighted by STD at each bin during
        %fitting process
        CCG_error = std(CCG_segment);
        fixed_param_flags = zeros(8,1); %by default, all 8 parameters will vary
        fixed_param_values = zeros(8,1); %override these values and flags to fix a parameter
        [pars{cond},freq(cond),Accounted_Var(cond)] = gaborgauss_sumfit([time', Raw_CCG(cond,:)'],raw,fixed_param_flags,fixed_param_values, CCG_error);
        fit_cross_correlogram(cond,:) = real (gaborgauss_sumfunc(time, pars{cond}) ) ;
        if pars{cond}(8)*2 <= 1/pars{cond}(5)*1000
            %gaussian contributing to the peak add gausian amplitude to
            %gabor amplitude, offset = baseline
            A(cond) = pars{cond}(2) + pars{cond}(6);
            O(cond) = pars{cond}(1);
        else
            %gaussian much broader than peak, add amplitude of gaussian to
            %baseline to get offset.
            A(cond) = pars{cond}(2);
            O(cond) = pars{cond}(1) + pars{cond}(6);
        end
        NC(cond) = A(cond)/O(cond) * 100;
        fit_gabor(cond,:) = pars{cond}(1) + pars{cond}(2)*exp( -(  abs(time - pars{cond}(3)) / pars{cond}(4)).^pars{cond}(7)).*cos(2*pi*pars{cond}(5)/1000 * (time - pars{cond}(3)) ); 
        fit_gaussian(cond,:) = pars{cond}(1) + pars{cond}(6) * exp( -((time - pars{cond}(3) )/ pars{cond}(8) ).^2);
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
    text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Peak t, Ht, A: %3d %6.3f %6.3f', cross_corr_peak_time(cond), cross_corr_peak_height(cond), cross_corr_peak_area(cond) );
    text(xpos,ypos, line,'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('NC: %6.3f, NCCs:  (%2d), %6.3f; (%2d), %6.3f; (%2d), %6.3f; (%2d), %6.3f; (%2d), %6.3f', NC(cond), peak_range, NCC(cond), peak_range2, NCC2(cond), peak_range3, NCC3(cond), peak_range4, NCC4(cond), peak_range5, NCC5(cond) );
    text(xpos,ypos, line, 'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Curve Fit Params: %6.1f %6.1f %6.3f %6.1f %6.1f %6.1f %6.3f %6.3f %6.1f', pars{cond} );
    text(xpos,ypos, line, 'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Perc. Var Accounted by Fit: %6.3f', Accounted_Var(cond) );
    text(xpos,ypos, line, 'FontSize',font_size);		ypos = ypos - bump_size;
    line = sprintf('Permutation NCC P value: %6.3f', permutation_NCC_p_values(cond) );
    text(xpos,ypos, line, 'FontSize',font_size);		ypos = ypos - bump_size;
    
    axis('off');
    hold on;
end	% of analyzing for each condition




%compare correlation strength using perm tet
if num_permutations2 > 0
    for cond = 1: length(comp_corr_conds)
        %normalize by number of segments in particular cond
        comp_corr_CCG{cond}  = comp_corr_CCG{cond} / length(mini_segment_start_times{comp_corr_conds(cond)});
        comp_corr_ACG1{cond}  = comp_corr_ACG1{cond}  / length(mini_segment_start_times{comp_corr_conds(cond)});
        comp_corr_ACG2{cond}  = comp_corr_ACG2{cond}  / length(mini_segment_start_times{comp_corr_conds(cond)});
        
        
    end
    
    %create vector that has identities for each trial over the two conditions
    %being compared
    trial_id = [1:num_reps(comp_corr_conds(1))  -(1:num_reps(comp_corr_conds(2))) ];
    
    for perm2 = 1:num_permutations2
        shuffle_index2 = trial_id(randperm(length(trial_id) ) );
        %choose correct number of trials for first cond
        trials1 = shuffle_index2(1:num_reps(comp_corr_conds(1)));
        %choose correct number for next cond
        trials2 = shuffle_index2(num_reps(comp_corr_conds(1)) + 1: end );
        
        perm_CCG1 = (sum(comp_corr_CCG{1}(trials1(trials1 > 0), : )) + sum(comp_corr_CCG{2}(-trials1(trials1 < 0), : ) ) )/ num_reps(comp_corr_conds(1));
        perm_ACG11 = (sum(comp_corr_ACG1{1} (trials1(trials1 > 0), : )) + sum(comp_corr_ACG1{2}(-trials1(trials1 < 0), : ) ) )/ num_reps(comp_corr_conds(1));
        perm_ACG21 = (sum(comp_corr_ACG2{1} (trials1(trials1 > 0), : )) + sum(comp_corr_ACG2{2}(-trials1(trials1 < 0), : ) ) )/ num_reps(comp_corr_conds(1));
        ACG11_area = sum( perm_ACG11( find(time == 0 ) - peak_range / neural_bin_width : find(time == 0 ) + peak_range / neural_bin_width )  );
        ACG21_area = sum( perm_ACG21( find(time == 0 ) - peak_range / neural_bin_width : find(time == 0 ) + peak_range / neural_bin_width )  );
        CCG1_area = sum( perm_CCG1( find(time == 0 ) - peak_range / neural_bin_width : find(time == 0 ) + peak_range / neural_bin_width )  );
        
        NCC_perm1 = CCG1_area / (ACG11_area * ACG21_area).^0.5;
        
        perm_CCG2  = (sum(comp_corr_CCG{1}(trials2(trials2 > 0), : )) + sum(comp_corr_CCG{2}(-trials2(trials2 < 0), : ) ))/ num_reps(comp_corr_conds(2)); 
        perm_ACG12  = (sum(comp_corr_ACG1{1}(trials2(trials2 > 0), : )) + sum(comp_corr_ACG1{2}(-trials2(trials2 < 0), : ) ))/ num_reps(comp_corr_conds(2)); 
        perm_ACG22  = (sum(comp_corr_ACG2{1}(trials2(trials2 > 0), : )) + sum(comp_corr_ACG2{2}(-trials2(trials2 < 0), : ) ))/ num_reps(comp_corr_conds(2)); 
        ACG12_area = sum( perm_ACG12( find(time == 0 ) - peak_range / neural_bin_width : find(time == 0 ) + peak_range / neural_bin_width )  );
        ACG22_area = sum( perm_ACG22( find(time == 0 ) - peak_range / neural_bin_width : find(time == 0 ) + peak_range / neural_bin_width )  );
        CCG2_area = sum( perm_CCG2( find(time == 0 ) - peak_range / neural_bin_width : find(time == 0 ) + peak_range / neural_bin_width )  );
        NCC_perm2 = CCG2_area / (ACG12_area * ACG22_area).^0.5;
        NCC_perm_diff(perm2) = NCC_perm1 - NCC_perm2;
    end
    
    NCC_diff = NCC(comp_corr_conds(1)) -  NCC(comp_corr_conds(2));
    NCC_diff_p_value(1) = sum(NCC_diff < NCC_perm_diff)/num_permutations2;
    
    if (length(comp_corr_conds ) == 3)
        trial_id = [1:num_reps(comp_corr_conds(1))  -(1:num_reps(comp_corr_conds(3))) ];
                
        for perm2 = 1:num_permutations2
            shuffle_index2 = trial_id(randperm(length(trial_id) ) );
            %choose correct number of trials for first cond
            trials1 = shuffle_index2(1:num_reps(comp_corr_conds(1)));
            %choose correct number for next cond
            trials2 = shuffle_index2(num_reps(comp_corr_conds(1)) + 1: end );
            
            perm_CCG1 = (sum(comp_corr_CCG{1}(trials1(trials1 > 0), : )) + sum(comp_corr_CCG{3}(-trials1(trials1 < 0), : ) ) )/ num_reps(comp_corr_conds(1));
            perm_ACG11 = (sum(comp_corr_ACG1{1} (trials1(trials1 > 0), : )) + sum(comp_corr_ACG1{3}(-trials1(trials1 < 0), : ) ) )/ num_reps(comp_corr_conds(1));
            perm_ACG21 = (sum(comp_corr_ACG2{1} (trials1(trials1 > 0), : )) + sum(comp_corr_ACG2{3}(-trials1(trials1 < 0), : ) ) )/ num_reps(comp_corr_conds(1));
            ACG11_area = sum( perm_ACG11( find(time == 0 ) - peak_range / neural_bin_width : find(time == 0 ) + peak_range / neural_bin_width )  );
            ACG21_area = sum( perm_ACG21( find(time == 0 ) - peak_range / neural_bin_width : find(time == 0 ) + peak_range / neural_bin_width )  );
            CCG1_area = sum( perm_CCG1( find(time == 0 ) - peak_range / neural_bin_width : find(time == 0 ) + peak_range / neural_bin_width )  );
            
            NCC_perm1 = CCG1_area / (ACG11_area * ACG21_area).^0.5;
            
            perm_CCG2  = (sum(comp_corr_CCG{1}(trials2(trials2 > 0), : )) + sum(comp_corr_CCG{3}(-trials2(trials2 < 0), : ) ))/ num_reps(comp_corr_conds(3)); 
            perm_ACG12  = (sum(comp_corr_ACG1{1}(trials2(trials2 > 0), : )) + sum(comp_corr_ACG1{3}(-trials2(trials2 < 0), : ) ))/ num_reps(comp_corr_conds(3)); 
            perm_ACG22  = (sum(comp_corr_ACG2{1}(trials2(trials2 > 0), : )) + sum(comp_corr_ACG2{3}(-trials2(trials2 < 0), : ) ))/ num_reps(comp_corr_conds(3)); 
            ACG12_area = sum( perm_ACG12( find(time == 0 ) - peak_range / neural_bin_width : find(time == 0 ) + peak_range / neural_bin_width )  );
            ACG22_area = sum( perm_ACG22( find(time == 0 ) - peak_range / neural_bin_width : find(time == 0 ) + peak_range / neural_bin_width )  );
            CCG2_area = sum( perm_CCG2( find(time == 0 ) - peak_range / neural_bin_width : find(time == 0 ) + peak_range / neural_bin_width )  );
            
            NCC_perm2 = CCG2_area / (ACG12_area * ACG22_area).^0.5;
            
            NCC_perm_diff(perm2) = NCC_perm1 - NCC_perm2;
        end
        NCC_diff = NCC(comp_corr_conds(1)) -  NCC(comp_corr_conds(3));
        NCC_diff_p_value(2) = sum(NCC_diff < NCC_perm_diff)/num_permutations2;
    end
    
        if (length(comp_corr_conds ) == 4)
        trial_id = [1:num_reps(comp_corr_conds(3))  -(1:num_reps(comp_corr_conds(4))) ];
        
        for perm2 = 1:num_permutations2
            shuffle_index2 = trial_id(randperm(length(trial_id) ) );
            %choose correct number of trials for first cond
            trials1 = shuffle_index2(1:num_reps(comp_corr_conds(1)));
            %choose correct number for next cond
            trials2 = shuffle_index2(num_reps(comp_corr_conds(1)) + 1: end );
            
            perm_CCG1 = (sum(comp_corr_CCG{3}(trials1(trials1 > 0), : )) + sum(comp_corr_CCG{4}(-trials1(trials1 < 0), : ) ) )/ num_reps(comp_corr_conds(3));
            perm_ACG11 = (sum(comp_corr_ACG1{3} (trials1(trials1 > 0), : )) + sum(comp_corr_ACG1{4}(-trials1(trials1 < 0), : ) ) )/ num_reps(comp_corr_conds(3));
            perm_ACG21 = (sum(comp_corr_ACG2{3} (trials1(trials1 > 0), : )) + sum(comp_corr_ACG2{4}(-trials1(trials1 < 0), : ) ) )/ num_reps(comp_corr_conds(3));
            ACG11_area = sum( perm_ACG11( find(time == 0 ) - peak_range / neural_bin_width : find(time == 0 ) + peak_range / neural_bin_width )  );
            ACG21_area = sum( perm_ACG21( find(time == 0 ) - peak_range / neural_bin_width : find(time == 0 ) + peak_range / neural_bin_width )  );
            CCG1_area = sum( perm_CCG1( find(time == 0 ) - peak_range / neural_bin_width : find(time == 0 ) + peak_range / neural_bin_width )  );
            
            NCC_perm1 = CCG1_area / (ACG11_area * ACG21_area).^0.5;
            
            perm_CCG2  = (sum(comp_corr_CCG{3}(trials2(trials2 > 0), : )) + sum(comp_corr_CCG{4}(-trials2(trials2 < 0), : ) ))/ num_reps(comp_corr_conds(4)); 
            perm_ACG12  = (sum(comp_corr_ACG1{3}(trials2(trials2 > 0), : )) + sum(comp_corr_ACG1{4}(-trials2(trials2 < 0), : ) ))/ num_reps(comp_corr_conds(4)); 
            perm_ACG22  = (sum(comp_corr_ACG2{3}(trials2(trials2 > 0), : )) + sum(comp_corr_ACG2{4}(-trials2(trials2 < 0), : ) ))/ num_reps(comp_corr_conds(4)); 
            ACG12_area = sum( perm_ACG12( find(time == 0 ) - peak_range / neural_bin_width : find(time == 0 ) + peak_range / neural_bin_width )  );
            ACG22_area = sum( perm_ACG22( find(time == 0 ) - peak_range / neural_bin_width : find(time == 0 ) + peak_range / neural_bin_width )  );
            CCG2_area = sum( perm_CCG2( find(time == 0 ) - peak_range / neural_bin_width : find(time == 0 ) + peak_range / neural_bin_width )  );
            
            NCC_perm2 = CCG2_area / (ACG12_area * ACG22_area).^0.5;
            
            NCC_perm_diff(perm2) = NCC_perm1 - NCC_perm2;
        end
        NCC_diff = NCC(comp_corr_conds(3)) -  NCC(comp_corr_conds(4));
        NCC_diff_p_value(2) = sum(NCC_diff < NCC_perm_diff)/num_permutations2;
    end

    
end



sub_max = 1.5 * max( max(CCG) );
sub_min = 1.5 * min( min(CCG) );

trial = 0;
% plotting loop

for cond =  1: num_conditions            
    if ~isempty(mini_segment_start_times{cond})
        subplot(num_rows, num_columns, cond*num_columns - 2 + plot_offset);
        phandle = plot(time', CCG(cond,:), 'r'  );
        XLim([-disp_range +disp_range]);
        YLim([sub_min sub_max]);
        if (cond == 1)
            title('Subtracted');
        end
        hold on         
                
        if (gabor_fit == 1)
            % plot gabor and gaussian parts of fit correlogram
            subplot(num_rows, num_columns, cond*num_columns - 1 + plot_offset);
            phandle = plot(time', Raw_CCG(cond,:) , 'r', time, fit_cross_correlogram(cond,:) , 'g',time', fit_gabor(cond,:) , 'b:', time', fit_gaussian(cond,:) , 'm:' );
            if (cond == 1)
                title('Fit Correlogram');
            end
        end
        
        hold on;
    end
    
    if plot_PSTH == 1    
        HIST_SCALE = 1.25 * max([max(max(hist_data2)) max(max(hist_data2B )) ]);
        subplot(num_rows, num_columns, cond*num_columns - 3 + plot_offset);
        phandle = plot(hist_time', hist_data2(cond,:), 'k-', hist_time', hist_data2B(cond,:), 'r-' );
        hold on;
        if (num_selected_segments(cond) > 1)
            xx = [mini_segment_start_times{cond}' mini_segment_start_times{cond}'];
            yy = [0 HIST_SCALE];
            if ~isempty(xx)
                plot(xx',yy', 'k--')
            end
            
            hold on;
            
            xx = [mini_segment_end_times{cond}' mini_segment_end_times{cond}'];
            yy = [0 max([max(hist_data2(cond,:)) max(hist_data2B(cond,:) ) ]) ];
            if ~isempty(xx)
                plot(xx',yy', 'g-.')
            end
        end
        
        hold on;
        if (cond == 1)
            title ('PSTH');
        end
        hold on; 
        XLim([hist_time(1) hist_time(end)]);
        YLim([0 HIST_SCALE]);
        
    end
       
end

%output cross correlation metrics
 



if (output == 1)
    eval(['save ' PATHOUT FILEOUT  '  CCG ACG1 ACG2 NCC_diff_p_value permutation_NCC5_p_values permutation_NCC4_p_values permutation_NCC3_p_values permutation_NCC2_p_values permutation_NCC_p_values NCC5 NCC4 NCC3 NCC2 peak_range2 peak_range3 peak_range4 peak_range5 Raw_CCG fit_cross_correlogram time fit_gaussian fit_gabor pars correlogram_norm_type minimum_segment_length num_selected_segments num_reps selected_spikes1 selected_spikes2 selected_FR1 selected_FR2 cross_corr_peak_time cross_corr_peak_height cross_corr_peak_area NCC NC num_permutations range peak_range peak_window spike_density_filter SpikeChanA SpikeChanB start_code stop_code BegTrial EndTrial StartOffset StopOffset StartEventBin StopEventBin PATH FILE'])
    
    %save figure
    FILEOUT = [FILE(1:i) 'fig'];        
    saveas (curr_fig, [PATHOUT FILEOUT]);
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function for computing correlogram matrices

function [correlogram, smooth_shuffle, corr_matrix, RawCCG] = compute_correlogram(spikesA, spikesB, time, neural_bin_width, norm_type)

num_reps = size(spikesA,1);
combin = 0;
triangle = ((size(spikesA, 2) )*neural_bin_width/1000) - abs(time)/1000;
corr_matrix = zeros(num_reps*num_reps, length(time) );
corr_matrix_GaborFit = zeros(num_reps*num_reps, length(time) );
smooth_shuffle = 0 * time;
correlogram = 0 * time;
RawCCG = 0 * time;
shuffle_GaborFit = 0 * time;
crosses = [];
crosses(linspace(1,(num_reps)^2, num_reps)) = 1;

switch norm_type
    case 0 %no normalization
        for trial1 = 1 : num_reps
            for trial2 = 1 : num_reps
                combin = combin + 1;
                corr_matrix(combin,:) = xcorr(spikesA(trial1, :), spikesB(trial2, : ), max(time) /neural_bin_width); 
                %convert to coincidences / sec by normalizing for finite trial length
                corr_matrix_GaborFit(combin,:) = corr_matrix(combin,:)./ triangle*size(spikesA,2)/1000;   
                corr_matrix(combin,:) = corr_matrix(combin,:)./ triangle;	
            end
        end  
        smooth_shuffle = mean(corr_matrix(~logical(crosses) , :)  );
        corr_matrix = corr_matrix - ones(size(corr_matrix,1), 1) * smooth_shuffle;
        correlogram = sum(corr_matrix(logical(crosses) , :)  );   
        RawCCG = sum(corr_matrix_GaborFit(logical(crosses),:)); 
        
    case 1 %normalize by geometric mean spike rate
        for trial1 = 1 : num_reps
            for trial2 = 1 : num_reps
                combin = combin + 1;
                corr_matrix(combin,:) = xcorr(spikesA(trial1, :), spikesB(trial2, : ), max(time) /neural_bin_width); 
                %convert to coincidences / sec by normalizing for finite trial length
                corr_matrix_GaborFit(combin,:) = corr_matrix(combin,:)./ triangle*size(spikesA,2)/1000;
                
                corr_matrix(combin,:) = corr_matrix(combin,:)./ triangle;	
                corr_matrix(combin,:) = corr_matrix(combin,:) /  ( ( (mean( spikesA(trial1, :) ) *1000  ) * (mean( spikesB(trial2, : ) )*1000 ) )^0.5 + 0.00001 );
            end
        end  
        % first obtain the general structure of the correlogram by creating smoothed shuffle
        % only from trial pairings that aren't cross correlograms
        smooth_shuffle = mean(corr_matrix(~logical(crosses) , :)  );
        corr_matrix = corr_matrix - ones(size(corr_matrix,1), 1) * smooth_shuffle;
        correlogram = mean(corr_matrix(logical(crosses) , :)  );   
        shuffle_GaborFit = mean(corr_matrix_GaborFit(~logical(crosses),:))*num_reps; 
        RawCCG = sum(corr_matrix_GaborFit(logical(crosses),:)) - shuffle_GaborFit + mean(shuffle_GaborFit); 
        
    case 2
        %normalize by minimum firing rate
        for trial1 = 1 : num_reps
            for trial2 = 1 : num_reps
                combin = combin + 1;
                corr_matrix(combin,:) = xcorr(spikesA(trial1, :), spikesB(trial2, : ), max(time) /neural_bin_width); 
                %convert to coincidences / sec by normalizing for finite trial length
                corr_matrix_GaborFit(combin,:) = corr_matrix(combin,:)./ triangle*size(spikesA,2)/1000;
                corr_matrix(combin,:) = corr_matrix(combin,:)./ triangle;	
                corr_matrix(combin,:) = corr_matrix(combin,:) / (min( [(mean( spikesA(trial1, :) ) *1000  ) (mean( spikesB(trial2, : ) )*1000 ) ] ) + 0.00001 );
            end
        end  
        % first obtain the general structure of the correlogram by creating smoothed shuffle
        % only from trial pairings that aren't cross correlograms
        smooth_shuffle = mean(corr_matrix(~logical(crosses) , :)  );
        corr_matrix = corr_matrix - ones(size(corr_matrix,1), 1) * smooth_shuffle;
        correlogram = mean(corr_matrix(logical(crosses) , :)  );   
        RawCCG = sum(corr_matrix_GaborFit(logical(crosses),:)); 
        
    case 3 %normalize by product of mean spike rate
        for trial1 = 1 : num_reps
            for trial2 = 1 : num_reps
                combin = combin + 1;
                corr_matrix(combin,:) = xcorr(spikesA(trial1, :), spikesB(trial2, : ), max(time) /neural_bin_width); 
                %convert to coincidences / sec by normalizing for finite trial length
                corr_matrix_GaborFit(combin,:) = corr_matrix(combin,:)./ triangle*size(spikesA,2)/1000;
                corr_matrix(combin,:) = corr_matrix(combin,:)./ triangle;	
                corr_matrix(combin,:) = corr_matrix(combin,:) / (mean( spikesA(trial1, :) ) *1000 * mean( spikesB(trial2, : ) )*1000  + 0.00001 );
            end
        end  
        % first obtain the general structure of the correlogram by creating smoothed shuffle
        % only from trial pairings that aren't cross correlograms
        smooth_shuffle = mean(corr_matrix(~logical(crosses) , :)  );
        corr_matrix = corr_matrix - ones(size(corr_matrix,1), 1) * smooth_shuffle;
        correlogram = mean(corr_matrix(logical(crosses) , :)  );   
        RawCCG = sum(corr_matrix_GaborFit(logical(crosses),:)); 
        
end      


