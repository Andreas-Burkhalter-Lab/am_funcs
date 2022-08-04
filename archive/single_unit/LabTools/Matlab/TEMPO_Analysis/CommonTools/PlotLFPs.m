function PlotLFPs(data, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol)

	% STILL NEED TO REVISE TO GET ORDINATE INTO SPIKES/S  - 9/25/00 - looks like scaling is good
	% check scaling for resample function
    line_types = ['b-'; 'r-'; 'g-'; 'k-'; 'g-'; 'g-'; 'g-'; 'g-'];
    Tempo_Defs;
    ProtocolDefs;
    
switch SpikeChan
    case 4 %first LFP channel
        neural_db = LFP_DB;
        SpikeChanA = 1;
    case 5 %second LFP channel
        neural_db = LFP_DB;
        SpikeChanA = 2;
    case 1
        neural_db = LFP_DB;
        SpikeChanA = 1;
    case 2
        neural_db = LFP_DB;
        SpikeChanA = 2;
end

switch SpikeChan2
    case 4 %first LFP channel
        neural_db = LFP_DB;
        SpikeChanB = 1;   
    case 5 %second LFP channel
        neural_db = LFP_DB;
        SpikeChanB = 2;
    case 1
        neural_db = LFP_DB;
        SpikeChanB = 1;
    case 2
        neural_db = LFP_DB;
        SpikeChanB = 2;
end

    low_band_low_cutoff = 0;
    low_band_high_cutoff = 8;

    band1_low_cutoff = 8;  %hz
    band1_high_cutoff = 30; %hz

    band2_low_cutoff = 20; %h
    band2_high_cutoff = 70; %Hz

    band3_low_cutoff = 70; %Hz
    band3_high_cutoff = 100; %Hz

    high_band_low_cutoff = 100;
    high_band_high_cutoff = 200;
   
    output = 1;
    
   %need uniform duration of analysis across ALL trials - use minimum difference between start and stop 
   %event bins and add offsets to find interval of analysis
   AnalInterval = min( (StopEventBin + StopOffset) - (StartEventBin + StartOffset) ) ;
   PreAnalInterval = 500 - mod(AnalInterval,500);
   PostAnalInterval = 500;
   
   %first, we'll need the bin_width (in sec) of our spike raster + events log
    h = data.htb_header{LFP_DB};	%for convenience
   LFP_bin_width = (h.skip + 1) / (h.speed_units / h.speed);
     
   h = data.htb_header{EVENT_DB};	%for convenience
   event_bin_width = (h.skip + 1) / (h.speed_units / h.speed);
       
   [conditions, unique_conds, param, num_conditions] = RegenerateConditionList(data, BegTrial, EndTrial, PATH, FILE, Protocol);

   num_modality = size(unique_conds,2) -1;
 
   %now, remove trials from hor_disp and spike_rates that do not fall between BegTrial and EndTrial
   trials = 1:size(conditions,2);												% a vector of trial indices
   select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );  
   
   start_eventbin = 1;   
   stop_eventbin = size(data.event_data, 2);   
   align_eventbin = find(data.event_data(1,:,1) == StartCode);    

   % convert start and stop times (for data range) to spike bins
   total_spike_bins = AnalInterval + PreAnalInterval;
   bin_width = 2;  % in ms
   num_reduced_bins = floor(total_spike_bins/bin_width);
   
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
        %had to introduce rounding here.  Sync pulses recorded every ms
        %while lfps are recorded every 2 ms.
        
        lfp = data.lfp_data(SpikeChanA, round ( ((StartEventBin(reps(1)  ) - PreAnalInterval): (LFP_bin_width*1000):StartEventBin( reps(1) ) + AnalInterval - 1 + PostAnalInterval)/(LFP_bin_width*1000) ) , reps(1)); 
        lfp2 = data.lfp_data(SpikeChanB, round ( ((StartEventBin(reps(1)  ) - PreAnalInterval): (LFP_bin_width*1000):StartEventBin( reps(1) ) + AnalInterval - 1 + PostAnalInterval)/(LFP_bin_width*1000) ) , reps(1)); 

        %     spike_rates(cond) = sum(data.lfp_data(SpikeChan, floor( StartEventBin(reps(1)  )/(LFP_bin_width * 1000) ) : floor ( ( StartEventBin( reps(1) ) + AnalInterval - 1/(AnalInterval*LFP_bin_width)  )  ), reps(1)  ) ); 
        rms(1, cond) = ( sum( (data.lfp_data(SpikeChanA, round( ((StartEventBin(reps(1)  ) - PreAnalInterval): (LFP_bin_width*1000): StartEventBin( reps(1) ) + AnalInterval - 1 + PostAnalInterval) / (LFP_bin_width * 1000)  ) , reps(1))) .^2 ) )^0.5; 
        rms(2,cond) =  ( sum( (data.lfp_data(SpikeChanB, round( ((StartEventBin(reps(1)  ) - PreAnalInterval): (LFP_bin_width*1000): StartEventBin( reps(1) ) + AnalInterval - 1 + PostAnalInterval) / (LFP_bin_width * 1000)  ) , reps(1))) .^2 )  )^0.5; 

%        full_rect_lfp(:,cond) = data.lfp_data(SpikeChanA, round( ((StartEventBin(reps(1)  ) - PreAnalInterval): (LFP_bin_width*1000): StartEventBin( reps(1) ) + AnalInterval - 1 + PostAnalInterval) / (LFP_bin_width * 1000)  ) , reps(1).^2); 

        
        hist_data = lfp;
        [freq, ampl] = FourierTransform_1D(linspace(0,LFP_bin_width*length(lfp),length(lfp)), lfp, length(lfp), 1,0);
        fourier_ampl(cond,:, 1) = ampl;
        [freq, ampl] = FourierTransform_1D(linspace(0,LFP_bin_width*length(lfp),length(lfp)), lfp2, length(lfp), 1,0);
        fourier_ampl(cond,:, 2) = ampl;

        low_band(cond,:, 1)  = sum( fourier_ampl(cond, (freq > low_band_low_cutoff) & (freq <= low_band_high_cutoff), 1 ) ); 
        band1(cond,:, 1)  = sum( fourier_ampl(cond, (freq > band1_low_cutoff) & (freq <= band1_high_cutoff), 1 ) ); 
        band2(cond,:, 1)  = sum( fourier_ampl(cond, (freq > band2_low_cutoff) & (freq <= band2_high_cutoff), 1 ) ); 
        band3(cond,:, 1)  = sum( fourier_ampl(cond,  (freq > band3_low_cutoff) & (freq <= band3_high_cutoff), 1 ) ); 
        high_band(cond,:, 1)  = sum(  fourier_ampl(cond, (freq > high_band_low_cutoff) & (freq <= high_band_high_cutoff), 1 ) ); 

        hist_dataB = lfp2;
        low_band(cond,:, 2)  = sum( fourier_ampl(cond,  (freq > low_band_low_cutoff) & (freq <= low_band_high_cutoff), 2 ) ); 
        band1(cond,:, 2)  = sum( fourier_ampl(cond, (freq > band1_low_cutoff) & (freq <= band1_high_cutoff), 2 ) ); 
        band2(cond,:, 2)  = sum( fourier_ampl(cond, (freq > band2_low_cutoff) & (freq <= band2_high_cutoff), 2 ) ); 
        band3(cond,:, 2)  = sum( fourier_ampl(cond,  (freq > band3_low_cutoff) & (freq <= band3_high_cutoff), 2 ) ); 
        high_band(cond,:, 2)  = sum(  fourier_ampl(cond,  (freq > high_band_low_cutoff) & (freq <= high_band_high_cutoff), 2 ) ); 
    
		for trial = 2: num_reps
      	    start_bin = StartEventBin(reps(trial)) + StartOffset;
            lfp = data.lfp_data(SpikeChanA, round( (start_bin - PreAnalInterval : (LFP_bin_width*1000): start_bin + AnalInterval - 1 + PostAnalInterval)/(LFP_bin_width * 1000)   ), reps(trial));   
            lfp2 = data.lfp_data(SpikeChanB, round( (start_bin - PreAnalInterval : (LFP_bin_width*1000): start_bin + AnalInterval - 1 + PostAnalInterval)/(LFP_bin_width * 1000)   ), reps(trial));   

            hist_data = hist_data + lfp;
           %         spike_rates(cond) = (spike_rates(cond)*(num_reps -1) + sum(data.lfp_data(SpikeChan, start_bin : start_bin + AnalInterval - 1, reps(trial)) )/(AnalInterval*LFP_bin_width) )/num_reps; 
            rms(1, cond) =  rms(1, cond) + (sum(data.lfp_data(SpikeChanA, round(  (start_bin - PreAnalInterval : (LFP_bin_width*1000): start_bin + AnalInterval - 1 + PostAnalInterval)/(LFP_bin_width * 1000)  ), reps(trial)).^2    ) )^0.5; 
            [freq, ampl] = FourierTransform_1D(linspace(0,LFP_bin_width*length(lfp),length(lfp)), lfp, length(lfp), 1,0);

            low_band(cond,:, 1)  = low_band(cond,:, 1) + sum( ampl( (freq > low_band_low_cutoff) & (freq <= low_band_high_cutoff) ) ); 
            band1(cond,:, 1)  = band1(cond,:, 1) + sum( ampl( (freq > band1_low_cutoff) & (freq <= band1_high_cutoff) ) ); 
            band2(cond,:, 1)  = band2(cond,:, 1) + sum( ampl( (freq > band2_low_cutoff) & (freq <= band2_high_cutoff) ) ); 
            band3(cond,:, 1)  = band3(cond,:, 1) + sum( ampl( (freq > band3_low_cutoff) & (freq <= band3_high_cutoff) ) ); 
            high_band(cond,:, 1)  = high_band(cond,:, 1) + sum( ampl( (freq > high_band_low_cutoff) & (freq <= high_band_high_cutoff) ) ); 
            fourier_ampl(cond,:, 1) = fourier_ampl(cond,:, 1)  + ampl;

            
            hist_dataB = hist_dataB + lfp2;
           %         spike_rates(cond) = (spike_rates(cond)*(num_reps -1) + sum(data.lfp_data(SpikeChan, start_bin : start_bin + AnalInterval - 1, reps(trial)) )/(AnalInterval*LFP_bin_width) )/num_reps; 
            rms(2, cond) =  rms(2, cond) +  ( sum(  data.lfp_data(SpikeChanB, round(  (start_bin - PreAnalInterval : (LFP_bin_width*1000): start_bin + AnalInterval - 1 + PostAnalInterval)/(LFP_bin_width * 1000)  ), reps(trial)).^2    )^0.5 ); 
            [freq, ampl] = FourierTransform_1D(linspace(0,LFP_bin_width*length(lfp),length(lfp)), lfp2, length(lfp), 1,0);
            low_band(cond,:, 2)  = low_band(cond,:, 2) + sum( ampl( (freq > low_band_low_cutoff) & (freq <= low_band_high_cutoff) ) ); 
            band1(cond,:, 2)  = band1(cond,:, 2) + sum( ampl( (freq > band1_low_cutoff) & (freq <= band1_high_cutoff) ) ); 
            band2(cond,:, 2)  = band2(cond,:, 2) + sum( ampl( (freq > band2_low_cutoff) & (freq <= band2_high_cutoff) ) ); 
            band3(cond,:, 2)  = band3(cond,:, 2) + sum( ampl( (freq > band3_low_cutoff) & (freq <= band3_high_cutoff) ) ); 
            high_band(cond,:, 2)  = high_band(cond,:, 2) + sum( ampl( (freq > high_band_low_cutoff) & (freq <= high_band_high_cutoff) ) ); 
            fourier_ampl(cond,:, 2) = fourier_ampl(cond,:, 2)  + ampl;

            
            
      end
      rms(:,cond) = rms(:,cond)/num_reps;
      fourier_ampl(cond,:, :) = fourier_ampl(cond,:, :)/num_reps;
      
      hist_data2(cond,:) = hist_data;   
      hist_data2B(cond,:) = hist_dataB;   

      scale(cond) = max(hist_data2(cond,:)  )  ;

      low_band(cond,:)  = low_band(cond,:)/num_reps;

      band1(cond,:)  = band1(cond,:)/num_reps; 
      band2(cond,:)  = band2(cond,:)/num_reps; 
      
      band3(cond,:)  = band3(cond,:)/num_reps;
      high_band(cond,:)  = high_band(cond,:)/num_reps;

      
   end   

    time = (-PreAnalInterval:(LFP_bin_width*1000):AnalInterval + PostAnalInterval - 1);

   num_columns = floor(num_conditions/5.01) + 1;
   if (num_columns > 3)
      num_columns = 3;
   end
   num_rows = floor(num_conditions/num_columns) + 1;
      
 %  	align_time = find(event_data{curr_condition, num_reps{curr_condition}} == align_code);
 %  marker_time = find(event_data{curr_condition, num_reps{curr_condition}} == marker_code);
   
     
   figure;
   set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [0 20 500 600], 'Name', 'PSTH Display');
   num_rows = num_conditions + 1;
    num_columns = 3;
    plot_offset = num_columns * 1;

   HIST_SCALE = 2.0;
   hist_max = max(hist_data2);	% a vector with max for each spike channel
   if ( max(hist_max) >= HIST_SCALE )		%if the histo for any spike is too large, rescale 'em all
      HIST_SCALE = round(max(hist_max));
   end

   trial = 0;
   % plotting loop
   for cond = 1: num_conditions 
        subplot(num_rows, num_columns, cond*num_columns - 2 + plot_offset);
        phandle = plot(time', hist_data2(cond,:), line_types(1), time', hist_data2B(cond,:), line_types(2));
   %    YLim([0 HIST_SCALE]);
        XLim([time(1) time(length(time))]);
        hold on  
        
        %set(phandle,'LineWidth', 1.5);
        xx = [0 0];
   		xx2 = [AnalInterval*LFP_bin_width AnalInterval*LFP_bin_width];
	    yy = [0 HIST_SCALE];
   		hold on;
   		plot(xx, yy, 'k-');
   		plot(xx2, yy,'k--');
		hold on;
        YLim( [min(min(hist_data2) ) max(max(hist_data2) ) ] );
        linetext =[];
         for modality = 1: num_modality
            linetext = [linetext param{modality}, num2str(unique_conds(cond,modality)), ' '];
         end   
         text('Position',[time(1) max(max(hist_data2) ) - 0.1],'String', linetext,'FontSize',8);         
       % gamma - 0-4 hz theta - 4-8 hz alpha 8-12 hz ? beta > 12 gamma ? ?
    
         subplot(num_rows, num_columns, cond*num_columns - 1 + plot_offset);
         phandle = plot(freq, fourier_ampl(cond,:, 1), line_types(1),  freq, fourier_ampl(cond,:, 2), line_types(2));       
         XLim([0 200]);        
            
   end
   
%     unique_dirs = unique_conds(2:end,2);
%     % plot tuning curve
%     subplot(ceil(num_rows/2),num_columns, num_columns);
%     plot(unique_dirs, band1(2:end) );
%     xlabel('Direction');
%     hold on;
%     null_y = [band1(1) band1(1)];
%     null_x = [min(unique_dirs) max(unique_dirs)];
% 	plot(null_x, null_y, 'k--');
%    
%    subplot(ceil(num_rows/2),num_columns, num_columns*2);
%     plot(unique_dirs, band2(2:end) );
%     xlabel('Direction');
%     null_y = [band2(1) band2(1)];
%     null_x = [min(unique_dirs) max(unique_dirs)];
% 	hold on;
% 	plot(null_x, null_y, 'k--');
% 
%    subplot(ceil(num_rows/2),num_columns, num_columns*3);
%     plot(unique_dirs, band3(2:end) );
%     xlabel('Direction');
%     null_y = [band3(1) band3(1)];
%     null_x = [min(unique_dirs) max(unique_dirs)];
% 	hold on;
% 	plot(null_x, null_y, 'k--');
% 
%    subplot(ceil(num_rows/2),num_columns, num_columns*4);
%     plot(unique_dirs, high_band(2:end) );
%     xlabel('Direction');
%     null_y = [high_band(1) high_band(1)];
%     null_x = [min(unique_dirs) max(unique_dirs)];
% 	hold on;
% 	plot(null_x, null_y, 'k--');
% 


    %output tuning curve metrics
    if (output == 1)
        i = size(PATH,2) - 1;
        while PATH(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
            i = i - 1;
        end   
        PATHOUT = [PATH(1:i) 'Analysis\LFP\'];
        i = size(FILE,2) - 1;
        while FILE(i) ~='.'
            i = i - 1;
        end
        FILEOUT = [FILE(1:i) 'lfp'];
        eval(['save ' PATHOUT FILEOUT  ' low_band band1 band2 band3 high_band rms   low_band_low_cutoff low_band_high_cutoff band1_low_cutoff band1_high_cutoff band2_low_cutoff band2_high_cutoff band3_low_cutoff band3_high_cutoff high_band_low_cutoff high_band_high_cutoff SpikeChanA SpikeChanB StartCode StopCode BegTrial EndTrial StartOffset StopOffset StartEventBin StopEventBin PATH FILE'])
    end

