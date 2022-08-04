%Edited on 1/20/08 by RLS to incorporate a seperate 'Pooling' function for use
%with calculating latency

function Pool_PSTH(bin_size, window_size, data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
line_types = ['b-'; 'r-'; 'g-'; 'k-'; 'g-'; 'g-'; 'g-'; 'g-'];
Tempo_Defs;
ProtocolDefs;

TEMPO_Defs;		%needed for defines like IN_T1_WIN_CD
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01
Path_Defs;

symbols = {'ko' 'k*' 'go' 'mo' 'b*' 'r*' 'g*' 'c*'};
lines = {'k-' 'm-' 'g-' 'b-' 'b--' 'r--' 'g--' 'c--'};
colors = {'k.' 'b.'};
global pool_data1 pool_data2 pool_size pool_trials 


%bin_matrix = [0:0.001:3]; %1 msec bins
%bin_centers = [0.0005:0.001:2.9995];

bin_matrix = [0:bin_size:3]; %2 msec bins
bin_centers = [bin_size/2:bin_size:3-(bin_size/2)];
bin_centers = [bin_centers, 3];


%figure;
%set(gcf,'PaperPosition', [.51 .51 20 27.5], 'Position', [450 50 500 573], 'Name', 'PSTH Graph');
%hold on
[hits1] = histc(pool_data1/1000, bin_matrix);
tally = [];
for i = 1:length(hits1);
    if i == 1;
        tally(i) = hits1(i);
    else
        tally(i) = hits1(i) + tally(i-1);
    end
end


[ave_psth] = hits1/pool_trials;
g = fspecial('gaussian', 3, 2);
gauss_filter = imfilter(ave_psth, g);


    
%recent_change = gauss_filter(1);
%no_change_steps = 0;
%for i = 2:length(gauss_filter);
%    if recent_change ~= gauss_filter(i);
%        if no_change_steps ~= 0;
%            slope(i) = (gauss_filter(i) - recent_change)/(.001*no_change_steps);
%            for j = 1:no_change_steps;
%                slope(i-j) = slope(i) - (j*(slope(i)- slope(i-no_change_steps-1))/no_change_steps);
%            end
%        else
%            slope(i) = (gauss_filter(i) - recent_change)/.001
%        end
%        recent_change = gauss_filter(i);
%        no_change_steps = 0;
%    else
%        no_change_steps = no_change_steps + 1;
%    end
%end
%diff_y = diff(tally);
%diff_x = diff(bin_centers);
%for i = 1:length(diff_x);
%    slope(i) = diff_y(i)/diff_x(i);
%end
%slope(length(slope)+1) = slope(length(slope));
%figure;
%set(gcf,'PaperPosition', [.51 .51 20 27.5], 'Position', [450 50 500 573], 'Name', 'PSTH');
%hold on
%plot(bin_centers, hits1, lines{1})
%figure;
%set(gcf,'PaperPosition', [.51 .51 20 27.5], 'Position', [450 50 500 573], 'Name', 'Average PSTH Graph');
%hold on
%plot(bin_centers, ave_hits1, lines{1})
figure;
set(gcf,'PaperPosition', [.51 .51 20 27.5], 'Position', [450 50 500 573], 'Name', 'Gauss Graph');
hold on
subplot(20, 1, 2:20);
plot(bin_centers, gauss_filter, lines{1})
subplot(20, 1, 1);
delay(window_size, bin_centers, gauss_filter, PATH, FILE);
%figure;
%set(gcf,'PaperPosition', [.51 .51 20 27.5], 'Position', [450 50 500 573], 'Name', 'Spikes v Time Graph');
%hold on
%plot(bin_centers, tally, lines{1})
%figure;
%set(gcf,'PaperPosition', [.51 .51 20 27.5], 'Position', [450 50 500 573], 'Name', 'Ave. PSTH');
%hold on
%plot(bin_centers, ave_psth, lines{1})

if ~isempty(pool_data2)
    [hits2] = histc(pool_data2/1000, bin_matrix);
    [ave_hits2] = hits2/5;
    [gauss_pool2] = (1/2*sqrt(2*pi))*exp(-ave_hits2/8);
    plot(bin_centers,gauss_pool2,lines{2})
end
title('Title')

return;

