%Edited on 06/19/06 by CMA to remove all references to SpikeChan2 or
%spikes2, which produce the error "??? Undefined function or variable
%'SpikeChan2'." when run with our data.  The original file is backed up in
%the folder Chris's Lab Tools Backup located on the desktop of this
%machine.

function SF_PSTH(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);

line_types = ['b-'; 'r-'; 'g-'; 'k-'; 'g-'; 'g-'; 'g-'; 'g-'];
Tempo_Defs;
ProtocolDefs;
data.neuron_params
pref_sf = data.neuron_params(PREFERRED_SPATIAL_FREQ, 1);
sf = data.gratings_params(GRAT_SPATIAL_FREQ,:,PATCH1);
diff_sf = pref_sf-sf;
diff_sf = abs(diff_sf);
diff_sf(diff_sf==min(diff_sf))=1;
diff_sf(diff_sf~=min(diff_sf))=0;
diff_sf(diff_sf~=0)=1;

total_spikes = [];
count=1;
figure
hold on
for trial=1:length(diff_sf);
    if diff_sf(trial)==1;
        count = count+1;
        spikes =  find(data.spike_data(SpikeChan,(StartEventBin(trial) + StartOffset - 500):(StopEventBin(trial) + StopOffset + 500),trial) == 1);
        plot (spikes/1000,count*ones(1,length(spikes)),'k.');
        total_spikes = cat(2,total_spikes,spikes);
        hold on;
    end
end
stim_bar = [1:3000];
on_off = ones(size(stim_bar));
on_off(1:500)=0;
on_off(2500:3000)=0;
plot(stim_bar/1000,on_off,'k')
axis([0 3 0 (count+1)])
xlabel('Time (s)')
ylabel('Trial')
title('Optimized Spatial Frequency Spike Rasters')

total_spikes = total_spikes/1000;
figure
hist(total_spikes, 60);
xlabel('Time (s)'), ylabel('Activity')
title('Grating Spatial Freqency PSTH')

%Save the total_spikes array to a cumulative data file so data can be
%analyzed for optimized SF accross multiple cells in each region of the
%brain.
BASE_PATH = ['C:\LabTools\Matlab\TEMPO_Analysis\'];
output = 0;
if output ==1;
    outfile = [BASE_PATH 'ProtocolSpecific\GratingSpatialFreq\ALSpikes_OptimizedSF.dat'];
    fid = fopen(outfile, 'a');
    dlmwrite(outfile, total_spikes, '-append')
    fclose(fid);
end
