%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% EventTimes.m - This function plots event times for each trial (in diff colors)
%	BJP - 4/26/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55

function EventTimes(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset,  StartEventBin, StopEventBin, PATH, FILE);

line_types = ['bo'; 'ro'; 'go'; 'co'; 'ko'; 'mo'; 'bd'; 'rd'; 'gd'; 'cd'; 'kd'; 'md'; 'b*'; 'r*'; 'g*'; 'c*'; 'k*'; 'm*'];
Tempo_Defs;
ProtocolDefs;



h = data.htb_header{EVENT_DB};	%for convenience
event_bin_width = (h.skip + 1) / (h.speed_units / h.speed);

figure;
set(gcf,'Name', 'Events Plot');


for trial = BegTrial:EndTrial
   event_times = find(data.event_data(1,:,trial) ~=0);   

   for marker = 1:length(event_times)
      phandle = plot(event_bin_width * event_times(marker), trial, line_types( data.event_data(1, event_times(marker), trial)  )   );
      hold on;
   end
   
   if size(data.spike_data, 1) > 1
   	%add synch pulse times for start and end frame
   	StimStart = min(find(data.spike_data(3,:,trial) == 1) );
   	StimStop = max(find(data.spike_data(3,:,trial) == 1) );
   	phandle = plot(event_bin_width * StimStart, trial, line_types( 1   ) );
   	hold on;
   	phandle = plot(event_bin_width * StimStop, trial, line_types( 2   ) );
		hold on;   
	end
      
end

title(['Events Plot:  Filename: ', PATH, FILE]);
XLABEL('Time (sec)');
YLABEL('Trial Number');
YLim([0 trial]);
XLim([0 event_bin_width * size(data.event_data,2)]);

%legend('TRIAL_START_CD', 'FP_ON_CD', 'IN_FIX_WIN_CD', 'VSTIM_ON_CD', 'VSTIM_OFF_CD','TARGS_ON_CD' 'SACCADE_BEGIN_CD','IN_T1_WIN_CD', 'IN_T2_WIN_CD', 'BROKE_FIX_CD', 'BROKE_VERG_CD', 'SUCCESS_CD', 'REWARD_CD', 'PUNISH_CD', 'TRIAL_END_CD',	'MICROSTIM_ON_CD', 'MICROSTIM_OFF_CD');


