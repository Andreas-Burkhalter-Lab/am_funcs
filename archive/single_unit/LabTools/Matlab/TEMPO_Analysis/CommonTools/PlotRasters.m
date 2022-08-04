%Edited on 06/19/06 by CMA to remove all references to SpikeChan2 or
%spikes2, which produce the error "??? Undefined function or variable
%'SpikeChan2'." when run with our data.  The original file is backed up in
%the folder Chris's Lab Tools Backup located on the desktop of this
%machine.

function PlotRasters(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);

line_types = ['b-'; 'r-'; 'g-'; 'k-'; 'g-'; 'g-'; 'g-'; 'g-'];
Tempo_Defs;
ProtocolDefs;

pref_sf = data.neuron_params(PREFERRED_SPATIAL_FREQ, 1);
sf = data.gratings_params(GRAT_SPATIAL_FREQ,:,PATCH1);
diff_sf = pref_sf-sf;
diff_sf = abs(diff_sf);
diff_sf(diff_sf==min(diff_sf))=1;
diff_sf(diff_sf~=min(diff_sf))=0;
diff_sf(diff_sf~=0)=1;

count=1;
figure
hold on
for trial=1:length(diff_sf);
    if diff_sf(trial)==1;
        count = count+1;
        spikes =  find(data.spike_data(SpikeChan,(StartEventBin(trial) + StartOffset - 500):(StopEventBin(trial) + StopOffset + 500),trial) == 1);
        plot (spikes/1000,count*ones(1,length(spikes)),'k.');
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

% figure
% hold on
% for trial=1:length(diff_sf);
%         count = count+1;
%         spikes =  find(data.spike_data(SpikeChan,:,trial) == 1);
%         plot (spikes,count*ones(1,length(spikes)),'k.');
%         hold on;
% end
% 