%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% ISIHist.m - This function plots first order interspike interval histograms
%	BJP - 4/30/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55

function ISIHist(data, Protocol, Analysis, SpikeChan, SpikeChan2, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset,  StartEventBin, StopEventBin, PATH, FILE);

line_types = ['bo'; 'ro'; 'go'; 'co'; 'ko'; 'mo'; 'bd'; 'rd'; 'gd'; 'cd'; 'kd'; 'md'; 'b*'; 'r*'; 'g*'; 'c*'; 'k*'; 'm*'];
Tempo_Defs;
ProtocolDefs;

h = data.htb_header{EVENT_DB};	%for convenience
event_bin_width = (h.skip + 1) / (h.speed_units / h.speed);
curr_fig = figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [0 20 500 600], 'Name', 'Interspike Interval Histograms');
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
 
%regenerate condition list
[conditions, unique_conds, param, num_conditions] = RegenerateConditionList(data, BegTrial, EndTrial, PATH, FILE, Protocol);  
num_modality = length(param);
    
num_rows = num_conditions + 2;
num_columns = 2;
plot_offset = num_columns * 2;
   
%now, remove trials from hor_disp and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:size(conditions,2);												% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );  
   
AnalInterval = min(StopEventBin - StartEventBin) + StartOffset + StopOffset;
spikes1 = zeros(1,AnalInterval);
spikes2 = zeros(1,AnalInterval);
   
for cond = 1: num_conditions
   
    SetTrials(1,:) = (conditions(1,:) == unique_conds(cond,1) & select_trials); 	    
    for modality = 1: num_modality           	
        NextSetTrials = (conditions(modality,:) == unique_conds(cond,modality) & select_trials); 
      	SetTrials = SetTrials & NextSetTrials;				
    end	 
     
    reps = find(SetTrials==1); 
    num_reps = length(reps);
     
    spikes1 = find (data.spike_data(SpikeChan, : ,reps(1) ) == 1);
    spikes2 = find (data.spike_data(SpikeChan, : ,reps(1) ) == 1);
    
	spikes1 = spikes1 (spikes1 >= StartEventBin(reps(1)) + StartOffset);
    spikes1 = spikes1 (spikes1 <= StartEventBin(reps(1)) + StartOffset + AnalInterval); 
	spikes2 = spikes2 (spikes2 >= StartEventBin(reps(1)) + StartOffset);
    spikes2 = spikes2 (spikes2 <= StartEventBin(reps(1)) + StartOffset + AnalInterval); 
    zero_bin = find (data.spike_data(SpikeChan, : ,reps(1) ) > 1);
    zero_bin2 = find (data.spike_data(SpikeChan2, : ,reps(1) ) > 1);
    isi{cond} = diff(spikes1);
    isi2{cond} = diff(spikes2);
    
    for trial = 2: num_reps
 	    spikes1 = find (data.spike_data(SpikeChan, : ,reps(trial) ) == 1);
 	    spikes2 = find (data.spike_data(SpikeChan2, : ,reps(trial) ) == 1);
      	start_bin = StartEventBin(reps(trial)) + StartOffset;
        spikes1 = spikes1 (spikes1 >= StartEventBin(reps(trial)) + StartOffset);
      	spikes1 = spikes1 (spikes1 <= StartEventBin(reps(trial)) + StartOffset + AnalInterval); 
        spikes2 = spikes2 (spikes2 >= StartEventBin(reps(trial)) + StartOffset);
      	spikes2 = spikes2 (spikes2 <= StartEventBin(reps(trial)) + StartOffset + AnalInterval); 
        zero_bin = [zero_bin find(data.spike_data(SpikeChan, : ,reps(1) ) > 1)];
        zero_bin2 = [zero_bin2 find(data.spike_data(SpikeChan2, : ,reps(1) ) > 1)];
      	isi{cond} = [isi{cond} diff(spikes1)]; 
      	isi2{cond} = [isi2{cond} diff(spikes2)]; 
        if ~isempty(zero_bin)
            isi{cond} = [isi{cond} zeros(1,length(zero_bin) )];
        end   
        if ~isempty(zero_bin2)
            isi2{cond} = [isi2{cond} zeros(1,length(zero_bin2) )];
        end   
        
    end  
end

for cond =  1: num_conditions            
    subplot(num_rows, num_columns, cond*num_columns - 1 + plot_offset);
    hist(isi{cond}, max(isi{cond}));
    linetext =[];   
    hold on       
    if cond == 1
        title(['First-order ISI Plots:  Filename: ', PATH, FILE]);
    end
    if cond == ceil(num_conditions/2)
        YLABEL('Count');
    end
    if cond == num_conditions
        XLABEL('Time (sec)');
    end
    xlim([0 max(isi{cond}) ]);
  
    subplot(num_rows, num_columns, cond*num_columns + plot_offset);
    hist(isi2{cond}, max(isi2{cond}));
    linetext =[];   
    hold on       
    if cond == 1
        title(['First-order ISI Plots:  Filename: ', PATH, FILE]);
    end
    if cond == ceil(num_conditions/2)
        YLABEL('Count');
    end
    if cond == num_conditions
        XLABEL('Time (sec)');
    end
    xlim([0 max(isi2{cond}) ]);

    
    
%    xlim([0 50 ]);
end

output = 1;

    %output tuning curve metrics
    if (output == 1)
        i = size(PATH,2) - 1;
        while PATH(i) ~='\'	%Analysis directory is one branch below Raw Data Dir
            i = i - 1;
        end   
        PATHOUT = [PATH(1:i) 'Analysis\ISI\'];
        i = size(FILE,2) - 1;
        while FILE(i) ~='.'
            i = i - 1;
        end
        FILEOUT = [FILE(1:i) 'isi'];
        eval(['save ' PATHOUT FILEOUT  ' isi isi2 BegTrial EndTrial StartOffset StopOffset StartEventBin StopEventBin PATH FILE'])
        FILEOUT = [FILE(1:i) 'fig'];        
        saveas (curr_fig, [PATHOUT FILEOUT]);
    end

