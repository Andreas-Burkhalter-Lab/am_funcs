function [binned_events, binned_spikes] = CEDread(PATH, FILE)
%global event_data spike_data binned_spikes 

event_string = '"CHANNEL"	"32"';
spike_string = '"CHANNEL"	"101"';
discriminator_string = '"CHANNEL"	"2"';

interval_before_trial = 1.0;
interval_after_trial = 3.0;

%read in the ced file
l = length(FILE);
if (FILE(l-3:l) == '.ced')	% .txt extension already there
	filename = [PATH FILE];   %the CED data file
elseif (FILE(l-3:l) == '.htb')	%.htb file extension, replace with .txt
	filename = [PATH FILE(1:l-4) '.ced'];   %the CED data file
end

file_events = [PATH FILE(1:l-4) '_events.txt'];
file_spikes = [PATH FILE(1:l-4) '_spikes.txt'];
file_dis_spikes = [PATH FILE(1:l-4) '_disc.txt'];

%files for input and output
fid = fopen(filename, 'rt');
fevents = fopen(file_events, 'wt');
fspikes = fopen(file_spikes, 'wt');
fdiscspikes = fopen(file_dis_spikes, 'wt');

%read in discriminator spikes
loopend = 0;
while loopend ~= 1
   in_line = fgetl(fid);
   test = strcmp(in_line, discriminator_string);
   if (test==1)
      loopend = 1;
      second_loop_end = 0;
      while second_loop_end ~= 1
         out_line = fgetl(fid);
         test2 = strcmp(out_line, event_string);
         if (test2 == 1)
            second_loop_end = 1;
         else
            fprintf(fevents, '%s \n', out_line);
      	end
   	end
   end
end
fclose(fdiscspikes);
disc_spike_times= textread(file_dis_spikes, '%f');
disc_spike_data = disc_spike_times;

loopend = 0;
while loopend ~= 1
   in_line = fgetl(fid);
   test = strcmp(in_line, event_string);
   if (test==1)
      loopend = 1;
      second_loop_end = 0;
      while second_loop_end ~= 1
         out_line = fgetl(fid);
         test2 = strcmp(out_line, spike_string);
         if (test2 == 1)
            second_loop_end = 1;
         else
            fprintf(fevents, '%s \n', out_line);
      	end
   	end
   end
end
fclose(fevents);
[event_times throw_away event_codes holder1 holder2 holder3] = textread(file_events, '%f %s %f %f %f %f');
event_data = [event_times event_codes];
 
%now the input file is at the correct place to read spike data
loopend = 0;
while loopend ~= 1
   out_line = fgetl(fid);
   if ~isstr(out_line)
      loopend = 1;
   else
      fprintf(fspikes, '%s \n', out_line);
   end
end
fclose(fspikes);
 
[spike_times throw_away spike_codes holder1 holder2 holder3] = textread(file_spikes, '%f %s %f %f %f %f');
spike_data = [spike_times spike_codes];

%now loop through and find a correct trial
%then find all the spikes within that trial
%all good trials will have a code 12 - successful trial
%then find the indices of the vis stim on of all successful trials
%then need to align trials to time code 4 - vis stim on
successful_trials = find(event_data(:, 2) == 12);
vis_stim_on_trials = zeros(length(successful_trials), 1);

j = 1;
for i = 1:length(event_data(:, 2))
    if (i==successful_trials(j))
        if (event_data(i-2, 2)==4)
            vis_stim_on_trials(j) = i-2;
            if (j < length(successful_trials))
                j = j + 1;
            end
        else
            disp('VIS_STIM_ON marker missing');
            break;
        end
    end
end

binned_spikes = zeros(length(successful_trials), 4000);
binned_events = zeros(length(successful_trials), 4000);
    

for i = 1:length(successful_trials)

    trial_begin = event_data(vis_stim_on_trials(i))-interval_before_trial;
    trial_end = event_data(vis_stim_on_trials(i))+interval_after_trial;  
    trial_length = trial_end-trial_begin;
    num_bins = round((trial_length/.001));
    
    single_trial_events = find((event_data(:, 1) >=trial_begin) & (event_data(:, 1) <=trial_end));
    event_times = event_data(single_trial_events) - trial_begin;
    event_index = floor(event_times * 1000);
    for j=1:length(single_trial_events)
        binned_events(i, event_index(j)) = event_data(single_trial_events(j), 2);
        binned_events(i, event_index(j));
    end
    %[val, temp_bins] = hist(event_data(single_trial_events, 2), num_bins)
    %binned_events(i, :) = temp_bins;
    
    single_trial = find((spike_data(:, 1) >= trial_begin) & (spike_data(:, 1) <= trial_end));
    temp_bins = hist(spike_data(single_trial), num_bins);
    binned_spikes(i, :) = temp_bins;
end
size(binned_spikes)

fclose(fid);

  