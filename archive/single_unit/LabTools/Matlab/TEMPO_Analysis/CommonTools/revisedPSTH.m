%*******************************************************************************************************************
%	PSTH.m - This function plots PSTHs for each spike channel (in diff. colors) for the current trial condition
%			All histos are rescaled together when one of them exceeds the limits of the Y axis.
%		BJP - 5/3/01
% AnalParams = column vector containing indices for analysis

function PSTH(data, SpikeChan, start_code, stop_code, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);

	% STILL NEED TO REVISE TO GET ORDINATE INTO SPIKES/S  - 9/25/00 - looks like scaling is good
	% check scaling for resample function

	line_types = ['b-'; 'r-'; 'g-'; 'k-'; 'g-'; 'g-'; 'g-'; 'g-'];
   Tempo_Defs;
   ProtocolDefs;
   
   %need uniform duration of analysis across ALL trials - use minimum difference between start and stop 
   %event bins and add offsets to find interval of analysis
   AnalInterval = min(StopEventBin - StartEventBin) + StartOffset + StopOffset;
   PreAnalInterval = 500 - mod(AnalInterval,500);
   PostAnalInterval = 500;
   
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
	            end     
   	   end 
      end   
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
   total_spike_bins = AnalInterval + PreAnalInterval;
   bin_width = 1;  % in ms
   num_reduced_bins = floor(total_spike_bins/bin_width);
	align_spikebins = floor (align_eventbin*(event_bin_width/spike_bin_width));	
   
   align_spikebins = 500; 
   
   
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
 		hist_data(1,:) = data.spike_data(SpikeChan, StartEventBin(reps(1)  ) - PreAnalInterval : StartEventBin( reps(1) ) + AnalInterval - 1 + PostAnalInterval, reps(1)); 
       
      spike_rates(cond) = sum(data.spike_data(SpikeChan, StartEventBin(reps(1)  ): StartEventBin( reps(1) ) + AnalInterval - 1, reps(1)) )/(AnalInterval*spike_bin_width); 
       

		for trial = 2: num_reps
      	 start_bin = StartEventBin(reps(trial)) + StartOffset;
          hist_data = hist_data +  data.spike_data(SpikeChan, start_bin - PreAnalInterval : start_bin + AnalInterval - 1 + PostAnalInterval, reps(trial));   
          
          spike_rates(cond) = (spike_rates(cond)*(num_reps -1) + sum(data.spike_data(SpikeChan, start_bin : start_bin + AnalInterval - 1, reps(trial)) )/(AnalInterval*spike_bin_width) )/num_reps; 
    
      end
      
      spike_rates(cond)
      
      % The next line requires the last two input arguments to be of interger proportions or else IT TAKES A LONG TIME
      %Resample of inputs ( 0 -> 1) gives outputs  ( 0 -> 0.025 ), with scaling factor of 40
      %resample of inputs (0 -> 6 ) gives outputs ( 0 -> .015) with scaling factor of 400
      % I think it is right now.  To check, I summed spikes between stim on and offset from hist_data.  This should equal the mean firing 
      % rate for each conditions (e.g. first condition = spont rate).
      % I compared this rate to the following:  I summed up all the values of hist_data2 between stim on and offset.  Then
      % I divided by the time (in seconds).  Since each bin in hist_data2 represented a four second window, I multiplied by 4
      % This gave a similar mean firing rate as for hist_data.
      
      hist_data2(cond,:) = 4*(resample(hist_data,num_reduced_bins, total_spike_bins)/num_reps);		
      %subsample the data after smoothing with a low-pass filter     
   	%each bin represents numbers of spikes in 4 sec window therefore, divide by four to get number of spikes/ms then multiply by 1000
      %to get spikes/s
%     hist_data2(cond,:) = 1000*hist_data2(cond,:)/4;
    fil = [.11 .11 .11 .11 .11 .11 .11 .11 .11];
    hist_data = conv(hist_data,fil);
    hist_data = hist_data( 5 :end - 4);
    hist_data2(cond,:) = hist_data;   
 % 5-23-01     scale(cond) = max(hist_data2(cond,floor(500/bin_width): floor(2000/bin_width)  )  );
      scale(cond) = max(hist_data2(cond,:)  )  ;
      
   end   

% 5-3-01 
%   time = [];
% 	for i=1:(num_reduced_bins)
% 		time(i) = (start_spikebin(1) + (total_spike_bins)*(i/num_reduced_bins) - align_spikebins)*spike_bin_width;
%   end   

time = (-PreAnalInterval:bin_width:AnalInterval + PostAnalInterval - 1)*spike_bin_width;

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
        
         %set(phandle,'LineWidth', 1.5);
         xx = [0 0];
   		xx2 = [AnalInterval*spike_bin_width AnalInterval*spike_bin_width];
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
 