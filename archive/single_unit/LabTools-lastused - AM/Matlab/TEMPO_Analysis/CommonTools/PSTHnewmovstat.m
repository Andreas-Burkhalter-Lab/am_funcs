%*******************************************************************************************************************
%	PSTH.m - This function plots PSTHs for each spike channel (in diff. colors) for the current trial condition
%			All histos are rescaled together when one of them exceeds the limits of the Y axis.
%			Greg DeAngelis, starting 9/1/99, edited from plothisto.m 
% AnalParams = column vector containing indices for analysis

function PSTH(data, SpikeChan, start_code, stop_code, BegTrial, EndTrial, StartOffset, StopOffset, StimStart, StimStop, PATH, FILE, Protocol);

	% STILL NEED TO REVISE TO GET ORDINATE INTO SPIKES/S  - 9/25/00 - looks like scaling is good
	% check scaling for resample function

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
   	   if size(munique(mod_values'),1) ~= 1
      	   num_modality = num_modality + 1;
         	conditions (num_modality,:) = mod_values; 
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
          end     %switch
   	   end % if loop
   end % for loop    
    % once have derived which conditions vary, get list of unique stimulus conditions, sort, and get total number
   unique_conds = munique(conditions');
   unique_conds = sortrows(unique_conds,1);
   num_conditions = size(unique_conds,1);
     
      
      
   else	% binding paradigms
      
      for modality = 1:1;
         for obj_num = 1:6;
 		  	   % get values for parameter 
		      mod_values = data.obj_params(modality,:,obj_num);
   		   if size(munique(mod_values'),1) ~= 1
  	    		   num_modality = num_modality + 1;
      	   	conditions (num_modality,:) = mod_values; 
       		   switch modality
    	 				case 1
      	   		param{num_modality} = ['Obj ' num2str(obj_num) ': '];    
	            end           
            end     
         end
      end
      
   conditions (end + 1,:) = data.misc_params(CONDITION, :);
   % once have derived which conditions vary, get list of unique stimulus conditions, sort, and get total number
   unique_conds = munique(conditions');
   unique_conds = sortrows(unique_conds,num_modality + 1)
   num_conditions = size(unique_conds,1);
      
      
   end
         
 
   %now, remove trials from hor_disp and spike_rates that do not fall between BegTrial and EndTrial
   trials = 1:size(conditions,2);												% a vector of trial indices
   select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );  
   
   start_eventbin = 1;   
   stop_eventbin = size(data.event_data, 2);   
   align_eventbin = find(data.event_data(1,:,1) == start_code);    

   % convert start and stop times (for data range) to spike bins
   start_spikebin = floor (start_eventbin*(event_bin_width/spike_bin_width));		
   stop_spikebin = floor (stop_eventbin*(event_bin_width/spike_bin_width));
   total_spike_bins = stop_spikebin(1) - start_spikebin(1) + 1;   % should be same for all trials
   bin_width = 20;  % in ms
   num_reduced_bins = floor(total_spike_bins/bin_width);
   align_spikebins = floor (align_eventbin*(event_bin_width/spike_bin_width));	  
   
   raw_hist = zeros(num_conditions,total_spike_bins);

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
          
          
 
          
          
      % accumulate data for all reps in given condition
     	for trial = 2: num_reps
         hist_data = hist_data +  data.spike_data(SpikeChan, start_spikebin : stop_spikebin, reps(trial));   
         %spikerate = 1000*sum(data.spike_data(SpikeChan, start_spikebin : stop_spikebin, reps(trial)))/(stop_spikebin - start_spikebin)
      end
      
      % average by the number of repetitions
      hist_data = hist_data/num_reps;
      
		% scaling = 50
      % The next line requires the last two input arguments to be of interger proportions or else IT TAKES A LONG TIME
      %Resample of inputs ( 0 -> 1) gives outputs  ( 0 -> 0.025 ), with scaling factor of 40
      %resample of inputs (0 -> 6 ) gives outputs ( 0 -> .015) with scaling factor of 400
      % I think it is right now.  To check, I summed spikes between stim on and offset from hist_data.  This should equal the mean firing 
      % rate for each conditions (e.g. first condition = spont rate).
      % I compared this rate to the following:  I summed up all the values of hist_data2 between stim on and offset.  Then
      % I divided by the time (in seconds).  Since each bin in hist_data2 represented a four second window, I multiplied by 4
      % This gave a similar mean firing rate as for hist_data.
      
 %      hist_data2(cond,:) = hist_data;
     
      hist_data2(cond,:) = 4*(resample(hist_data,num_reduced_bins, total_spike_bins)/num_reps);		
      %subsample the data after smoothing with a low-pass filter     
   	%each bin represents numbers of spikes in 4 sec window therefore, divide by four to get number of spikes/ms then multiply by 1000
      %to get spikes/s
      hist_data2(cond,:) = 1000*hist_data2(cond,:)/4;
         
      scale(cond) = max(hist_data2(cond,floor(find(data.event_data(1,:,1) == start_code)/bin_width): floor(find(data.event_data(1,:,1) == stop_code)/bin_width)));
      

      %spike_rate = computed average spike rate using resampled data
      spike_rate(cond) = sum(hist_data2(cond,50:125))/75;

      
   end   
     
   %%% Code for plotting composite PSTH
  
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

	%calculate max spike rate from moving trials at pref disp       
%   moving_max = scale(and(unique_conds(:,2) == pref_hdisp, unique_conds(:,1) >0));
 
   %% grab moving PSTH for preferred disparity only
%   moving_PSTH = hist_data2(and(unique_conds(:,2) == pref_hdisp, unique_conds(:,1) >0),:)/moving_max;
    moving_PSTH = hist_data2(logical(unique_conds(:,1) >0),:);
    moving_PSTH = sum(moving_PSTH)/size(moving_PSTH,1); 

   %% grab stationary PSTH for preferred disparity from only
%   stationary_PSTH = hist_data2(and(unique_conds(:,2) == pref_hdisp, unique_conds(:,1) == 0),:)/moving_max;
   stationary_PSTH = hist_data2(logical(unique_conds(:,1) == 0),:);
   stationary_PSTH = sum(stationary_PSTH)/size(stationary_PSTH,1);
   
   %% grab spontaneous PSTH
   spontaneous_PSTH = hist_data2(unique_conds(:,1) == data.one_time_params(NULL_VALUE),:);

   %normalize by peak value within VSTIM period
   stationary_max = max(stationary_PSTH( floor(StimStart(1)/bin_width) : floor(StimStop(1)/bin_width )  ) );
   
   moving_max = max(moving_PSTH( floor(StimStart(1)/bin_width) : floor(StimStop(1)/bin_width )  )  );    
   
   if stationary_max > moving_max
       moving_PSTH = moving_PSTH/stationary_max;
       stationary_PSTH = stationary_PSTH/stationary_max;
       spontaneous_PSTH = spontaneous_PSTH/stationary_max;
   else 
       moving_PSTH = moving_PSTH/moving_max;
       stationary_PSTH = stationary_PSTH/moving_max;
       spontaneous_PSTH = spontaneous_PSTH/moving_max;
   end
   
   time = [];
 	for i=1:(num_reduced_bins)
 		time(i) = (start_spikebin(1) + (total_spike_bins)*(i/num_reduced_bins) - align_spikebins)*spike_bin_width;
   end   
    
   num_columns = floor(num_conditions/5.01) + 1;
   if (num_columns > 3)
      num_columns = 3;
   end
   num_rows = floor(num_conditions/num_columns) + 1;
      
 %  	align_time = find(event_data{curr_condition, num_reps{curr_condition}} == align_code);
 %  marker_time = find(event_data{curr_condition, num_reps{curr_condition}} == marker_code);
   
     
   figure;
   set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [0 20 700 800], 'Name', 'PSTH Display');
   
   HIST_SCALE = 2.0;
   hist_max = max(hist_data2);	% a vector with max for each spike channel
   if ( max(hist_max) >= HIST_SCALE )		%if the histo for any spike is too large, rescale 'em all
      HIST_SCALE = round(max(hist_max));
   end

   trial = 0;
   % plotting loop
   for i = 1: num_columns 
      for j = 1: num_rows
         trial = trial + 1;
         if trial <= num_conditions
         subplot(num_rows, num_columns, trial);
         phandle = plot(time', hist_data2(trial,:), line_types(1));
         YLim([0 HIST_SCALE]);
         XLim([time(1) time(length(time))]);
         linetext =[];
         for modality = 1: num_modality
            linetext = [linetext param{modality}, num2str(unique_conds(trial,modality)), ' '];
         end   
         text('Position',[time(1) (HIST_SCALE - .4)],'String', linetext,'FontSize',8);         
         hold on  
         start_eventbin = find(data.event_data(1,:,trial) == start_code);   
  			stop_eventbin = find(data.event_data(1,:,trial) == stop_code);   
         % convert to spike bins
      	start_time = (-align_spikebins + StartOffset + floor (start_eventbin*(event_bin_width/spike_bin_width)))/1000;		
      	stop_time = (-align_spikebins + StopOffset + floor (stop_eventbin*(event_bin_width/spike_bin_width)))/1000;

         %set(phandle,'LineWidth', 1.5);
         xx = [start_time start_time];
   		xx2 = [stop_time stop_time];
		   yy = [0 HIST_SCALE];
   		hold on;
   		plot(xx, yy, 'k-');
   		plot(xx2, yy,'k--');
			hold on;
       
         end
     	end
   end
   
   subplot(num_rows,num_columns,1);
	text('Position',[time(1) (HIST_SCALE + HIST_SCALE/8)],'String', ['Filename: ', PATH, FILE],'FontSize',12);         
 
   composite = 1;
   
   if composite == 1
       
      figure
      plot(time, spontaneous_PSTH, time, stationary_PSTH, time, moving_PSTH);
      axis([min(time) max(time) 0 1]);
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
         stat_PSTH_out = zeros(1,length(stationary_PSTH));
         mov_PSTH_out = zeros(1,length(moving_PSTH));
         spont_PSTH_out = zeros(1,length(spontaneous_PSTH));
         pref_disp = [];
         site = [];
         time_out = [];
      end   
     
        if  ((length(time) < length(time_out) ) | (num_files == 0)  )
            stat_PSTH_out = (stat_PSTH_out(: , 1:length(time)  )*num_files + stationary_PSTH)/(num_files + 1);      
            mov_PSTH_out = (mov_PSTH_out(: , 1:length(time) ) *num_files + moving_PSTH)/ (num_files + 1);
            spont_PSTH_out = (spont_PSTH_out(: , 1:length(time) )*num_files + spontaneous_PSTH)/(num_files + 1);
            time_out = time;   
        elseif (length(time) > length(time_out))
            stat_PSTH_out = (stat_PSTH_out*num_files + stationary_PSTH(: , 1:length(time_out) ) )/(num_files + 1);      
            mov_PSTH_out = (mov_PSTH_out*num_files + moving_PSTH(: , 1:length(time_out) )  )/ (num_files + 1);
            spont_PSTH_out = (spont_PSTH_out*num_files + spontaneous_PSTH(: , 1:length(time_out) ) )/(num_files + 1);
            time_out = time(1:length(time_out) );   
        elseif (length(time) < length(time_out))
            stat_PSTH_out = (stat_PSTH_out*num_files + stationary_PSTH(: , 1:length(time) ) )/(num_files + 1);      
            mov_PSTH_out = (mov_PSTH_out*num_files + moving_PSTH(: , 1:length(time) )  )/ (num_files + 1);
            spont_PSTH_out = (spont_PSTH_out*num_files + spontaneous_PSTH(: , 1:length(time) ) )/(num_files + 1);
            time_out = time(1:length(time));   
        end     
           num_files = num_files + 1;
      
      pref_disp = [pref_disp; pref_hdisp]; 
            
      save (fileid, 'stat_PSTH_out', 'mov_PSTH_out', 'num_files', 'pref_disp', 'spont_PSTH_out', 'time_out');
    
%      figure
%      plot(time, spontaneous_PSTH, time, stationary_PSTH, time, moving_PSTH);
%      axis([min(time) max(time) 0 1]);
 %     title('Composite PSTH');
end % composite