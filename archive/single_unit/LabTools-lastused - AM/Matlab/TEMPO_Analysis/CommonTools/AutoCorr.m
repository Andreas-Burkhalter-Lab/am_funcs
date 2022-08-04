%*******************************************************************************************************************
% AutoCorr - Generates Auto Correlograms BJP 4/26/01
% 			Prints Autocorrelogram, shuffled correlogram, and auto - shuffled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function AutoCorr(data, Protocol, Analysis, SpikeChan, start_code, stop_code, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE);

    line_types = ['b-'; 'r-'; 'g-'; 'k-'; 'g-'; 'g-'; 'g-'; 'g-'];
   TEMPO_Defs;
   ProtocolDefs;

   neural_db = SPIKE_DB;
 %  neural_db = LFP_DB;
   
   %first, we'll need the bin_width (in sec) of our spike raster + events log
   h = data.htb_header{neural_db};	%for convenience
   neural_bin_width = (h.skip + 1) / (h.speed_units / h.speed) * 1000;  %bin width of neural signal in ms

   %need uniform duration of analysis across ALL trials - use minimum difference between start and stop 
   %event bins and add offsets to find interval of analysis
   AnalInterval = min(StopEventBin - StartEventBin) + StartOffset + StopOffset;
     
   num_reduced_bins = 100;
   align_time = 0; 

   range = 100;		
   time = -range:+range;
   %triangle correction factor for finite trial length
   triangle = (size(data.spike_data,2)*neural_bin_width/1000) - abs(time)/1000;
  
   trial = 0;
    
   [conditions, unique_conds, param, num_conditions] = RegenerateConditionList(data, BegTrial, EndTrial, PATH, FILE, Protocol);
   num_modality = length(param);
   
   %now, remove trials from hor_disp and spike_rates that do not fall between BegTrial and EndTrial
   trials = 1:size(conditions,2);												% a vector of trial indices
   select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );  
   
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
      histo(cond, :) = zeros(1, size(time,2));
      histo2(cond,:) = zeros(1, size(time,2));
      
      spikes = [];
      
      % Sequentially compute AutoCorrelograms averaged across trials in condition
	   for trial = 1:num_reps
   			start_spikebin = StartEventBin(reps(trial)) + StartOffset;
      	   switch neural_db
           case SPIKE_DB
                spikes(trial,:) = data.spike_data(SpikeChan, start_spikebin: start_spikebin + AnalInterval, reps(trial)); 
       	        histo(cond, :) = histo(cond,:) + xcorr(spikes(trial,:), spikes(trial,:), range)./ triangle;;		     
           case LFP_DB
                spikes(trial,:) = data.lfp_data(SpikeChan, floor(start_spikebin/neural_bin_width): floor( (start_spikebin)/neural_bin_width) + floor( (AnalInterval)/neural_bin_width )  , reps(trial)); 
        	    histo(cond, :) = histo(cond,:) + xcorr(spikes(trial,:), spikes(trial,:), range)./ triangle;;		      
           end     
            %histo(cond,find(time == 0)) = 0;         
            
      end  
         
      %shuffle trials      
      shuffle1 = randperm(num_reps);
      shuffle2 = randperm(num_reps);
      for trial = 1:num_reps
         %need this logic so that autocorrelograms are not computed during shuffling
      	if (shuffle1(trial) == shuffle2(trial)) 
         	others = [shuffle2(1:trial - 1) shuffle2(trial + 1:end)];
            shuffle2(trial) =  others(floor((num_reps - 1) * rand(1) ) + 1 ); 
         end
         histo2(cond, :) = histo2(cond,:) + xcorr(spikes(shuffle1(trial), :), spikes(shuffle2(trial), : ), range)./ triangle;		     
     	end    
            
	histo(cond, :) = histo(cond, :)/num_reps;
   histo2(cond, :) = histo2(cond, :)/num_reps;
   histo3(cond, :) = histo(cond,:) - histo2(cond,:);
	end	% of analyzing for each condition
   
   
   num_columns = floor(num_conditions/5.01) + 1;
   if (num_columns > 3)
      num_columns = 3;
   end
   num_rows = floor(num_conditions/num_columns) + 1;
   
   %num_rows = 4;
   %num_cols = 2;
   figure;
   set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [0 20 500 600], 'Name', 'Auto Correlograms');
   
   hist_max = max(histo);	% a vector with max for each spike channel
   auto_scale = max(hist_max);
   auto_min = min(min(histo));
   
   %% matrix containing composite corr data
   comp_hist = zeros(num_modality, size(histo,2));
  	%text('Position',[time(1) auto_scale + auto_scale/8],'String', ['Filename: ', PATH, FILE],'FontSize',12);         
	title(['Filename: ', PATH, FILE]);

   trial = 0;
   % plotting loop
   for i = 1: num_columns 
      for j = 1: num_rows
         trial = trial + 1;
	         if trial <= num_conditions
   	      	subplot(num_rows, num_columns, trial);
   	      	phandle = plot(time', histo(trial,:), line_types(1));
      	   	YLim([auto_min auto_scale]);
         		XLim([-range +range]);
         		linetext =[];
    		     	for modality = 1: num_modality 
         	   	linetext = [linetext param{modality}, num2str(unique_conds(trial,modality)), ' '];
       		  	end   
               title(linetext);   
               hold on  
            end   				                          
     	end
   end
     
   % plot shuffles
   num_columns = floor(num_conditions/5.01) + 1;
   if (num_columns > 3)
      num_columns = 3;
   end
   num_rows = floor(num_conditions/num_columns) + 1;
   
   %num_rows = 4;
   %num_cols = 2;
   figure;
   set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [0 20 500 600], 'Name', 'Shuffled Correlograms');
   
   hist_max = max(histo2);	% a vector with max for each spike channel
   auto_scale = max(hist_max);
   auto_min = min(min(histo2));
   
  	%text('Position',[time(1) auto_scale + auto_scale/8],'String', ['Filename: ', PATH, FILE],'FontSize',12);         
	title(['Filename: ', PATH, FILE]);

   subplot(num_rows,num_columns,1);
   trial = 0;
   % plotting loop
   for i = 1: num_columns 
      for j = 1: num_rows
         trial = trial + 1;
	         if trial <= num_conditions
   	      	subplot(num_rows, num_columns, trial);
   	      	phandle = plot(time', histo2(trial,:), line_types(1));
      	   	YLim([auto_min auto_scale]);
         		XLim([-range +range]);
         		linetext =[];
    		     	for modality = 1: num_modality 
         	   	linetext = [linetext param{modality}, num2str(unique_conds(trial,modality)), ' '];
       		  	end   
               title(linetext);   
               hold on  
            end  
     	end
   end
     
   subplot(num_rows,num_columns,1);
   text('Position',[time(1) auto_scale + auto_scale/8],'String', ['Filename: ', PATH, FILE],'FontSize',12);         


   % plot difference
   num_columns = floor(num_conditions/5.01) + 1;
   if (num_columns > 3)
      num_columns = 3;
   end
   num_rows = floor(num_conditions/num_columns) + 1;
   
  % num_rows = 4;
   %num_cols = 2;
   figure;
   set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [0 20 500 600], 'Name', 'Auto Correlograms (shuffled subtracted)');
   
   hist_max = max(histo3);	% a vector with max for each spike channel
   auto_scale = max(hist_max);
   auto_min = min(min(histo3));
   
  	%text('Position',[time(1) auto_scale + auto_scale/8],'String', ['Filename: ', PATH, FILE],'FontSize',12);         
	title(['Filename: ', PATH, FILE]);

   subplot(num_rows,num_columns,1);
   trial = 0;
   % plotting loop
   for i = 1: num_columns 
      for j = 1: num_rows
         trial = trial + 1;
	         if trial <= num_conditions
   	      	subplot(num_rows, num_columns, trial);
   	      	phandle = plot(time', histo3(trial,:), line_types(1));
      	   	YLim([auto_min auto_scale]);
         		XLim([-range +range]);
         		linetext =[];
    		     	for modality = 1: num_modality 
         	   	linetext = [linetext param{modality}, num2str(unique_conds(trial,modality)), ' '];
       		  	end   
               title(linetext);   
               hold on  
            end  
     	end
   end
     
   subplot(num_rows,num_columns,1);
   text('Position',[time(1) auto_scale + auto_scale/8],'String', ['Filename: ', PATH, FILE],'FontSize',12);         


   