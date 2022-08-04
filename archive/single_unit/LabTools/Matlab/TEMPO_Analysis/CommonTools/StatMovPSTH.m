%*******************************************************************************************************************
%	StatMovPSTH.m - This function compiles PSTH matrices for comparison of moving and
%           stationary trials.
%			BJP 7/16/01
% AnalParams = column vector containing indices for analysis

function StatMovPSTH(data, SpikeChan, start_code, stop_code, BegTrial, EndTrial, StartOffset, StopOffset, StimStart, StimStop, PATH, FILE, Protocol);

line_types = ['b-'; 'r-'; 'g-'; 'k-'; 'g-'; 'g-'; 'g-'; 'g-'];
Tempo_Defs;
ProtocolDefs;

%first, we'll need the bin_width (in sec) of our spike raster + events log
h = data.htb_header{SPIKE_DB};	%for convenience
spike_bin_width = (h.skip + 1) / (h.speed_units / h.speed);
     
h = data.htb_header{EVENT_DB};	%for convenience
event_bin_width = (h.skip + 1) / (h.speed_units / h.speed);
       
num_modality = 0;
% get unique conditions  
if (Protocol < 200)
    for modality = 1:17;
   	   % get values for parameter 
	    mod_values = data.dots_params(modality,:,PATCH1);
        mod_values = mod_values(~isnan(mod_values))';
   	    if ((~isempty(mod_values) ) & (size(munique(mod_values),1) ~= 1))
      	    num_modality = num_modality + 1;
            conditions (num_modality,:) = mod_values'; 
           	switch modality
                case 1
          	   		param{num_modality} = 'Dir: ';
      		    case 2                
         		    param{num_modality} = 'Speed: ';
                case 3
      	   	    	param{num_modality} = 'Coh: ';
   	   	        case 4
      	   	    	param{num_modality} = 'H Disp: ';
      		    case 5
         		    param{num_modality} = 'V Disp: ';
      		    case 6
         	    	param{num_modality} = 'Bin Corr: ';
      		    case 7
         	    	param{num_modality} = 'H Grad Mag: ';
      		    case 8
                 	param{num_modality} = 'H Grad Ang: ';     
                case 17
                    param{num_modality} = 'Patch Size: ';         
	        end     
         end 
    end   

    % once have derived which conditions vary, get list of unique stimulus conditions, sort, and get total number
    unique_conds = munique(conditions');
    unique_conds = sortrows(unique_conds,1);
    num_conditions = size(unique_conds,1);
         
end
         
 
%now, remove trials from hor_disp and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:size(conditions,2);												% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );  
   
start_eventbin = 1;   
stop_eventbin = size(data.event_data, 2);   

% convert start and stop times (for data range) to spike bins
start_spikebin = floor (start_eventbin*(event_bin_width/spike_bin_width));		
stop_spikebin = floor (stop_eventbin*(event_bin_width/spike_bin_width));
   
start_time = find(data.event_data(1,:,1) == start_code);    
stop_time = find(data.event_data(1,:,1) == stop_code);    
 
start_offset = 100; % ms prior to stimulus start
start_spikebin = start_time - start_offset;
stop_spikebin = stop_time;
time = (start_spikebin:1:stop_spikebin) - start_time + 1;   

spike_rate = zeros(1,num_conditions);
%analyze by unique condition
for cond = 1: num_conditions
    % These next few lines check for trials belonging to a single condition
    SetTrials(1,:) = (conditions(1,:) == unique_conds(cond,1) & select_trials); 	    
    for modality = 1: num_modality           	
        NextSetTrials = (conditions(modality,:) == unique_conds(cond,modality) & select_trials); 
      	SetTrials = SetTrials & NextSetTrials;				
    end	 
     
    reps = find(SetTrials==1); 
    num_reps = length(reps);
		   
	% initialize for first trial
 	hist_data(1,:) = data.spike_data(SpikeChan, start_spikebin : stop_spikebin, reps(1)); 
    spike_rate(cond) = spike_rate(cond) + data.spike_rates(reps(1) ) / num_reps;
     
 
    % accumulate data for all reps in given condition
    for trial = 2: num_reps
         hist_data = hist_data +  data.spike_data(SpikeChan, start_spikebin : stop_spikebin, reps(trial));   
         spike_rate(cond) = spike_rate(cond) + data.spike_rates(reps(trial) ) / num_reps;
    end
      
    % average by the number of repetitions
    hist_data2(cond,:) = hist_data/num_reps;
end   
       
    %% choose conditions that are not controls
    select_conds = logical(unique_conds(:,2) ~=LEYE_CONTROL & unique_conds(:,2) ~=REYE_CONTROL & unique_conds(:,2) ~=UNCORR_CONTROL & unique_conds(:,2) ~= data.one_time_params(NULL_VALUE));
  
    %% get spike rates for moving trials only
   moving_spike_rates = spike_rate(and(unique_conds(:,1) >0, select_conds));
   %% get spike rates for stationary trials only
   stationary_spike_rates = spike_rate(and(unique_conds(:,1) == 0, select_conds));

    %initialize marker for preferred disp
   pref_hdisp = -1;
   % choose preferred disparity by finding max spike rate for moving
   pref_hdisp = unique_conds(find(spike_rate == max(moving_spike_rates)),2)
 
%% grab moving PSTH for preferred disparity only
   moving_PSTH = hist_data2(and(unique_conds(:,2) == pref_hdisp, unique_conds(:,1) >0),:);

 
%% grab stationary PSTH for preferred disparity from only
stationary_PSTH = hist_data2(and(unique_conds(:,2) == pref_hdisp, unique_conds(:,1) == 0),:);

%% grab moving PSTH for preferred disparity only
%moving_PSTH = hist_data2(logical(unique_conds(:,1) >0),:);
%moving_PSTH = sum(moving_PSTH)/size(moving_PSTH,1); 

%% grab stationary PSTH for preferred disparity from only
%stationary_PSTH = hist_data2(logical(unique_conds(:,1) == 0),:);
%stationary_PSTH = sum(stationary_PSTH)/size(stationary_PSTH,1);
%% grab spontaneous PSTH
spontaneous_PSTH = hist_data2(unique_conds(:,1) == data.one_time_params(NULL_VALUE),:);

%normalize by peak value within VSTIM period
stationary_max = max(stationary_PSTH   );
   
moving_max = max(moving_PSTH );    
   
if stationary_max > moving_max
    moving_PSTH = moving_PSTH/stationary_max;
    stationary_PSTH = stationary_PSTH/stationary_max;
    spontaneous_PSTH = spontaneous_PSTH/stationary_max;
else 
    moving_PSTH = moving_PSTH/moving_max;
    stationary_PSTH = stationary_PSTH/moving_max;
    spontaneous_PSTH = spontaneous_PSTH/moving_max;
end
 
i = size(PATH,2) - 1;
while PATH(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
    i = i - 1;
end  
i = i - 1; 	%jump over slash
while PATH(i) ~='\'	%Go to Tempo directory
    i = i - 1;
end 
      
PATHOUT = [PATH(1:i)];
      
      
      
FILEOUT = 'PSTH.mat';  
fileid = [PATHOUT FILEOUT];
   
if exist(fileid, 'file') ~= 0
    load (fileid)
else 
    num_files = 0;
    stat_PSTH_out = []; %zeros(1,length(stationary_PSTH));
    mov_PSTH_out = []; %zeros(1,length(moving_PSTH));
    spont_PSTH_out = []; %zeros(1,length(spontaneous_PSTH));
    sites = [];
    time_out = time;
end   
     
if (length(time) >= length(time_out) )
    % new data longer or same legth than existing 
    stat_PSTH_out = [stat_PSTH_out;  stationary_PSTH(1:length(time_out) ) ];      
    mov_PSTH_out = [mov_PSTH_out; moving_PSTH(1:length(time_out) )];
    spont_PSTH_out = [spont_PSTH_out; spontaneous_PSTH(1:length(time_out) ) ];
    time_out = time_out;
else  %(max(time) < max(time_out) )
    % new data shorter than existing 
    stat_PSTH_out = [stat_PSTH_out(:,1:length(time)) ;  stationary_PSTH];      
    mov_PSTH_out = [mov_PSTH_out(:,1:length(time))  ; moving_PSTH];
    spont_PSTH_out = [spont_PSTH_out(:,1:length(time)) ; spontaneous_PSTH];
    time_out = time;
end

num_files = num_files + 1;
sites{num_files, 1} = FILE;      
save (fileid, 'stat_PSTH_out', 'mov_PSTH_out', 'num_files', 'spont_PSTH_out', 'sites', 'time_out');
      
figure
plot(time_out, mean(spont_PSTH_out,1), time_out, mean(stat_PSTH_out, 1), time_out, mean(mov_PSTH_out,1));
title('Average PSTH');