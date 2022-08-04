%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% JointPSTH.m - This function plots spike rasters for each spike channel and for each trial (in diff colors)
%	BJP - 4/30/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55

function JointPSTH(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset,  StimStart, StimStop, PATH, FILE);

line_types = ['bo'; 'ro'; 'go'; 'co'; 'ko'; 'mo'; 'bd'; 'rd'; 'gd'; 'cd'; 'kd'; 'md'; 'b*'; 'r*'; 'g*'; 'c*'; 'k*'; 'm*'];
Tempo_Defs;
ProtocolDefs;

if (SpikeChan == 1)
   SpikeChan1 = 1;
   SpikeChan2 = 2;
else
   SpikeChan1 = 2;
   SpikeChan2 = 1;
end

h = data.htb_header{EVENT_DB};	%for convenience
event_bin_width = (h.skip + 1) / (h.speed_units / h.speed);

figure;

  TrialLength = min(StimStop - StimStart);
   
   spikes1 = zeros(1,StimStart + TrialLength);
   spikes2 = spikes1;
   
for trial = BegTrial:EndTrial
   spikes1 = find (data.spike_data(SpikeChan1, StimStart(trial):StimStart(trial) + TrialLength ,trial) >= 1);
   
   spikes2 = find (data.spike_data(SpikeChan2,StimStart(trial):StimStart(trial) + TrialLength,trial) >= 1); 
    
   %spikes2 = spikes2(spikes2 > StimStart(trial) & spikes2 < StimStart(trial) + TrialLength)
   phandle = plot(event_bin_width * spikes1, trial*ones(1,length(spikes1)), 'b.' );
   hold on;
   phandle = plot(event_bin_width * spikes2, (trial + 0.5)*ones(1,length(spikes2)), 'r.' );
   
   
   
end

title(['Events Plot:  Filename: ', PATH, FILE]);
XLABEL('Time (sec)');
YLABEL('Trial Number');
YLim([0 trial]);
XLim([0 event_bin_width * size(data.event_data,2)]);

%legend('TRIAL_START_CD', 'FP_ON_CD', 'IN_FIX_WIN_CD', 'VSTIM_ON_CD', 'VSTIM_OFF_CD','TARGS_ON_CD' 'SACCADE_BEGIN_CD','IN_T1_WIN_CD', 'IN_T2_WIN_CD', 'BROKE_FIX_CD', 'BROKE_VERG_CD', 'SUCCESS_CD', 'REWARD_CD', 'PUNISH_CD', 'TRIAL_END_CD',	'MICROSTIM_ON_CD', 'MICROSTIM_OFF_CD');


