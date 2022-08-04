function delay(window_size, bin_centers, signal, PATH, FILE);
line_types = ['b-'; 'r-'; 'g-'; 'k-'; 'g-'; 'g-'; 'g-'; 'g-'];
Tempo_Defs;
ProtocolDefs;

TEMPO_Defs;		%needed for defines like IN_T1_WIN_CD
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01
Path_Defs;

symbols = {'ko' 'k*' 'go' 'mo' 'b*' 'r*' 'g*' 'c*'};
lines = {'k-' 'm-' 'g-' 'b-' 'b--' 'r--' 'g--' 'c--'};
colors = {'k.' 'b.'};

count = .5/bin_centers(1)/2;
pre_stim = [signal(1:count)];
post_stim = [signal(length(signal)-count:length(signal))];
unstim = [pre_stim, post_stim];
unstim_mean = mean(unstim);

%moving window
%window_size = .022;
flag = 0;
latency = 0;
for i = 1:length(signal);
    if i+((window_size/bin_centers(1)/2)-1) > length(signal);
        break
    end
    window = [signal(i:(i+((window_size/bin_centers(1)/2)-1)))];
    window_mean = mean(window);
    window_std = std(window);
    if window_mean > unstim_mean + 2*window_std;
        flag = flag + 1;
    else
        flag = 0;
    end
    if flag == 10;
        latency = bin_centers(i);
        break
    end
end


        
        
%i = 1;

%while bin_centers(i) <= .5 - bin_centers(1);
%    pre_stim = pre_stim + signal(i);
%    ave_pre_stim = pre_stim/i;
%    i = i + 1;
%end

%filtered_signal = signal - ave_pre_stim;

%[peak, peak_index] = max(filtered_signal);
%peak_time = bin_centers(peak_index);

%j = peak_index;

%for k = peak_index:-1:1
%    if (filtered_signal(k) >= .05 * peak) && (bin_centers(k) > .5);
%        latency = bin_centers(k);
%    end
%end

latency_time = latency - .5;

axis([0 100 0 100]);
axis('off');
font_size = 8;
bump_size = 60;
          
% type out stats onto graph
xpos = -10;   
ypos = 10;
temp = strcat(PATH, FILE(1:8));
temp(temp == '\') = '/';
% this prevents a stupid error from appearing on the screen
line = sprintf('File: %s', temp);
%Next line is the source of stupid error
text(xpos,ypos, line,'FontSize',font_size);
xpos = xpos + bump_size;
line = sprintf('Latency time (seconds): %f', latency_time);
text(xpos,ypos,line,'FontSize',font_size);	

output = 1;
if (output == 1)

    %------------------------------------------------------------------------
    %write out all relevant parameters to a cumulative text file, GCD 8/08/01
    %write out one line for each stimu_type for each neuron.
    outfile = [BASE_PATH 'Latency.txt'];
    printflag = 0;
    if (exist(outfile, 'file') == 0)    %file does not yet exist
        printflag = 1;
    end
    fid = fopen(outfile, 'a');
    if (printflag)
        fprintf(fid, 'Neuron\t\t Window Size\t Bin Size\t Latency\t');
        fprintf(fid, '\r\n');
        printflag = 0;
    end
    buff = sprintf('%s\t %5.3f\t\t %5.3f\t\t %5.5f\t', ...
        FILE(1:8), window_size, 2*bin_centers(1), latency_time);
    fprintf(fid, '%s', buff);
    fprintf(fid, '\r\n');
    fclose(fid);
    %------------------------------------------------------------------------
end
