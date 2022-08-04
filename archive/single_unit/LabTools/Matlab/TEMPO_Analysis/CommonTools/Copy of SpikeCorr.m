function [data, choptime] = SpikeCorr(data, windowflag, SpikeChan, SpikeChan2)

num_trials = size(data.spike_data);
num_trials = num_trials(3);

correl_total = 0;
for i=1:num_trials
   spike_data1 = data.spike_data(SpikeChan, :, i);
   spike_data2 = data.spike_data(SpikeChan2, :, i);
   [correl, lags] = xcorr(spike_data1, spike_data2, 50);
   correl_total = correl+correl_total;
end

[val, ind] = max(correl_total);
choptime = lags(ind);

%figure
%plot(lags, correl_total);

%spike_data1 = data.spike_data(SpikeChan, :);
%spike_data2 = data.spike_data(SpikeChan2, :);

%[correl, lags] = xcorr(spike_data1, spike_data2, 50);
%figure
%plot(lags, correl);

%[val, ind] = max(correl_total);
%choptime_short = lags_total(ind);

if windowflag == 0
    spike_data1 = data.spike_data(SpikeChan, :, :);
    spike_data2 = data.spike_data(SpikeChan2, :, :);

    sizem = size(spike_data1);
    num_trials = sizem(3);
    %load up all the trials from the spike_data
    for i = 1:num_trials
        spikes{i} = find(spike_data1(1, :, i) ~= 0);
        spikes2{i}= find(spike_data2(1, :, i) ~= 0);
    end

    for i = 1:num_trials %for each trial
        for j=1:length(spikes{i}) %for number of spikes in each trial
            %for each spike in the window discriminator
            if spikes{i}(j) - choptime > 0
                if spike_data2(1,spikes{i}(j) - choptime,i) ~= 0
                    spike_data2(1,spikes{i}(j) - choptime,i) = spike_data2(1,spikes{i}(j) - choptime,i)-1;
                end
            end   
    %        if spikes{i}(j) - choptime+1 > 0
    %           if spike_data2(1,spikes{i}(j) - choptime+1,i) ~= 0
    %                spike_data2(1,spikes{i}(j) - choptime+1,i) = spike_data2(1,spikes{i}(j) - choptime+1,i)-1;
    %            end
    %        end
    %        if spikes{i}(j) - choptime-1 > 0
    %            if spike_data2(1, spikes{i}(j) - choptime-1, i) ~= 0
    %                spike_data2(1, spikes{i}(j) - choptime-1, i) = spike_data2(1, spikes{i}(j) - choptime-1, i)-1;
    %            end
    %        end
        end
    end

    %temp1 = spike_data1(1,:);
    %temp2 = spike_data2(1,:);
    
    %[correl, lags] = xcorr(temp1, temp2, 50);
    %hold on
    %plot(lags, correl, 'r');
    
    data.spike_data(SpikeChan2, :, :) = spike_data2;
end