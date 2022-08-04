%Edited on 02/24/07 by CMA to perform fourier analysis of cell responses at
%optimal spatial frequency

function Fourier_PSTH(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, StartEventBin, StopEventBin, PATH, FILE, Protocol);
line_types = ['b-'; 'r-'; 'g-'; 'k-'; 'g-'; 'g-'; 'g-'; 'g-'];
Tempo_Defs;
ProtocolDefs;

TEMPO_Defs;		%needed for defines like IN_T1_WIN_CD
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01
Path_Defs;

symbols = {'ko' 'k*' 'go' 'mo' 'b*' 'r*' 'g*' 'c*'};
lines = {'k-' 'k--' 'g-' 'm-' 'b--' 'r--' 'g--' 'c--'};
colors = {'k.' 'b.'};

%get the column of values of SFs in gratings_params matrix
sf = data.gratings_params(GRAT_SPATIAL_FREQ,:,PATCH1);

%get indices of any NULL conditions (for measuring spontaneous activity)
null_trials = logical( (sf == data.one_time_params(NULL_VALUE)) );

%get the column of stimulus type values
stim_type = data.gratings_params(GRAT_TYPE,:,PATCH1);
unique_stim_type = munique(stim_type(~null_trials)');

%now, get the firing rates for all the trials
spike_rates = data.spike_rates(SpikeChan, :);

unique_sf = munique(sf(~null_trials)');

%now, remove trials that do not fall between BegTrial and EndTrial
trials = 1:length(sf);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );


total_spikes1 = [];
total_spikes2 = [];
count=1;
% figure(2);
% set(gcf,'PaperPosition', [.51 .51 20 27.5], 'Position', [450 50 500 573], 'Name', 'Optimized Spatial Frequency Raster');
dt = 0.001; %sampling at 1000 Hz, so the timestep is 0.001

%GAUSSFUNC Used by GAUSSFIT.
%  GAUSSFUNC assumes a function of the form
%
%	  y = q(1) + q(2) * exp(-0.5*((x - q(3))/q(4)).^2 )
%
%	thus q(1) is base rate, q(2) is amplitude, q(3) is center, and q(4) is size
sigma = 0.02;
x = [-4*sigma:dt:4*sigma];
q = [0, 1, 0, sigma];
z = q(1) + q(2) * exp(-0.5*((x - q(3))/ q(4)).^2);
% disp('length(z) =')
% length(z)

Y1_total = zeros(1,2000);
Y2_total = zeros(1,2000);
count1 = 0;
count2 = 0;

for j = 1:length(unique_stim_type);

    grating_select = logical( (stim_type == unique_stim_type(j)) );

    plot_x=sf(~null_trials & select_trials & grating_select);
    plot_y=spike_rates(~null_trials & select_trials & grating_select);

    %NOTE:  Inputs to PlotTuningCurve.m must be column vectors, not row vectors
    %because of munique().
    [px, py, perr] = PlotTuningCurve(plot_x', plot_y', symbols{j}, lines{j}, 1, 0);
    %NOTE: last argument=0 means just get output, no plot

    %py is the response and px is the stimulus value.  We are interested in
    %the cell's response at optimized conditions, therefore, we find the
    %max of py, and select the corresponding px as stimulus value that
    %results in optimized response.
    [max_response, index] = max(py);
    optimized_stimulus = px(index);

    totalmagnitude = zeros(1,2000);
    hold on
    for trial=1:length(sf);
        if (sf(trial) == optimized_stimulus & stim_type(trial) == unique_stim_type(j));
            count = count+1;

            %spikes =  data.spike_data(SpikeChan,:,trial);
            % from PlotSFRasters protocol
            spikes =  data.spike_data(SpikeChan,(StartEventBin(trial) + StartOffset):(StartEventBin(trial) + StartOffset + 1999),trial);
            spikes(spikes~=1) = 0;

            %filtered_spikes = conv(spikes,z);

            %Y = fft(filtered_spikes, 2000);
            %After convolution with the gaussian, the filtered spikes
            %vector will be longer than the original data vector, thus take
            %an FFT of the first 2000 points.
            %CMA 02/16/07

            %Used this to confirm that the 13th data point represents the
            %magnitude at frequency = 2 Hz
            %disp('f(5)')
            %f(5)

            if j == 1;
                if isempty(total_spikes1)
                    total_spikes1 = spikes;
                else
                    total_spikes1 = total_spikes1 + spikes;
                end
                %Y1_total = Y1_total + Y;
                %count1 = count1 + 1;
            elseif j == 2;
                if isempty(total_spikes2)
                    total_spikes2 = spikes;
                else
                    total_spikes2 = total_spikes2 + spikes;
                end
                %Y2_total = Y2_total + Y;
                %count2 = count2 + 2;
            else
                'Error: More than two grating types'
            end

        end
    end
end


filtered_spikes1 = conv(total_spikes1,z);
Y1 = fft(filtered_spikes1, 2000);
if ~isempty(total_spikes2)
    filtered_spikes2 = conv(total_spikes2,z);
    Y2 = fft(filtered_spikes2, 2000);
end

Y1_avg = Y1/max(Y1);
magnitude1 = abs(Y1_avg);

if ~isempty(total_spikes2)
    Y2_avg = Y2/max(Y2);
    magnitude2 = abs(Y2_avg);
end


% %Uncomment these lines to see a plot of the filtered spikes
% figure
% plot(filtered_spikes1)

f = (0:(2000-1))/2000/dt;


% %Plot the time domain representation of the data points with this:
% t = [1:1:2000]/1000;
% figure
% plot(t,filtered_spikes(1:length(t)))
% xlabel('Time (s)'), ylabel('Amplitude'), title('With Filtering')

% %Used this to confirm that the second data point and final data point are
% %equal
% magnitude(2)
% magnitude(length(magnitude))

figure
hold on
plot(0,0,'k')
plot(0,0,'k--')
plot(0,0,'ko')
legend('Grating 1', 'Grating 2', '0 Hz and 2 Hz Frequencies')
cla
plot(f(1:50),magnitude1(1:50),'k')
if ~isempty(total_spikes2)
    plot(f(1:50),magnitude2(1:50),'k--')
end

xlabel('Frequency (Hz)'), ylabel('Average Amplitude'), title(FILE), xlim([-1 25])
hold on
plot(f(5),magnitude1(5),'ko')
if ~isempty(total_spikes2)
    plot(f(5),magnitude2(5),'ko')
end
plot(f(1),magnitude1(1),'ko')
if ~isempty(total_spikes2)
    plot(f(1),magnitude2(1),'ko')
end

if ~isempty(total_spikes2)
    zero_Hz_frequency = [magnitude1(1), magnitude2(1)]./2;
    two_Hz_frequency = [magnitude1(5), magnitude2(5)];
    ratio = two_Hz_frequency./zero_Hz_frequency;
else
    zero_Hz_frequency = [magnitude1(1)]./2;
    two_Hz_frequency = [magnitude1(5)];
    ratio = two_Hz_frequency/zero_Hz_frequency;
end

output = 1;
if (output == 1)

    %------------------------------------------------------------------------
    %write out all relevant parameters to a cumulative text file, GCD 8/08/01
    %write out one line for each stimu_type for each neuron.
    outfile = [BASE_PATH 'ProtocolSpecific\GratingSpatialFreq\GratingSF_Fourier_Summary_AL.dat'];
    printflag = 0;
    if (exist(outfile, 'file') == 0)    %file does not yet exist
        printflag = 1;
    end
    fid = fopen(outfile, 'a');
    if (printflag)
        fprintf(fid, 'FILE\t\t\t PrfOR\t PrfSF\t PrfTF\t RFXctr\t RFYctr\t RFDiam\t GrType\t 0Hz_freq\t 2Hz_freq\t 2Hz/0Hz\t');
        fprintf(fid, '\r\n');
        printflag = 0;
    end
    for j = 1:length(unique_stim_type)
        buff = sprintf('%s\t %5.1f\t %5.2f\t %5.3f\t %5.2f\t %5.2f\t %5.2f\t %5.1f\t %5.4f\t %5.4f\t %5.4f\t', ...
            FILE, data.neuron_params(PREFERRED_ORIENTATION, 1), data.neuron_params(PREFERRED_SPATIAL_FREQ, 1), data.neuron_params(PREFERRED_TEMPORAL_FREQ, 1), data.neuron_params(RF_XCTR, 1), data.neuron_params(RF_YCTR, 1), data.neuron_params(RF_DIAMETER, 1),...
            unique_stim_type(j), zero_Hz_frequency(j), two_Hz_frequency(j), ratio(j));
        fprintf(fid, '%s', buff);
        fprintf(fid, '\r\n');
    end
    fclose(fid);
    %------------------------------------------------------------------------

end


