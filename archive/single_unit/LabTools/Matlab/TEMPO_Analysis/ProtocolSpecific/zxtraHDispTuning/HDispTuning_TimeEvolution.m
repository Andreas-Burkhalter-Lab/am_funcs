%-----------------------------------------------------------------------------------------------------------------------
%-- HDispTuning_TimeEvolution.m -- Plots horizontal disparity tuning as a function fo time
%--	GCD, starting 9/11/01
%-----------------------------------------------------------------------------------------------------------------------
function HDispTuning_TimeEvolution(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01
Path_Defs;

symbols = {'ko' 'k*' 'go' 'mo' 'b*' 'r*' 'g*' 'c*'};
lines = {'k-' 'k--' 'g-' 'm-' 'b--' 'r--' 'g--' 'c--'};

%get the column of values of horiz. disparity in the dots_params matrix
hor_disp = data.dots_params(DOTS_HDISP,:,PATCH1);

%get indices of any NULL conditions (for measuring spontaneous activity)
null_trials = logical( (hor_disp == data.one_time_params(NULL_VALUE)) );

%get the column of speed values
speed = data.dots_params(DOTS_SPEED,:,PATCH1);
unique_speed = munique(speed(~null_trials)');

%now, get the firing rates for all the trials 
spike_rates = data.spike_rates(SpikeChan, :);

%get indices of monoc. and uncorrelated controls
control_trials = logical( (hor_disp == LEYE_CONTROL) | (hor_disp == REYE_CONTROL) | (hor_disp == UNCORR_CONTROL) );

%now, remove trials from hor_disp and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(hor_disp);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

% Calculate spontaneous rates before looping through so can calculate DTI
null_rate = mean(data.spike_rates(SpikeChan, null_trials & select_trials));

figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [450 50 500 573], 'Name', 'Horizontal Disparity Tuning Curve');

window_hwid = 30; % half-width of moving window, ms
window_step = 20;  % increment of window position, ms

DDI = [];  var_term = []; DTI = []; PSTH_max = [];
for i=1:length(unique_speed)	%for each different speed value, plot a separate disparity tuning curve
    
    DT_data = [];
    time_slice = [];
    index = 1;
    for j = 0:window_step:1020 
        
        time_slice(index) = j;
        
        speed_select = logical( (speed == unique_speed(i)) );
        
        spike_rates = ComputeSpikeRates(data, length(hor_disp), StartCode, StartCode, j-window_hwid, j+window_hwid);
        % use the below line for a cumulative time plot
        %spike_rates = ComputeSpikeRates(data, length(hor_disp), StartCode, StartCode, 0, j);
        
        plot_x = hor_disp(speed_select & ~null_trials & ~control_trials & select_trials);
        plot_y = spike_rates(speed_select & ~null_trials & ~control_trials & select_trials); 
        
        %compute the Disparity Discrimination Index
        [DDI(i,index), var_term(i,index)] = Compute_DDI(plot_x, plot_y);
        
        %NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
        hold on;
        [px, py, perr, pmax, pmin] = PlotTuningCurve(plot_x', plot_y', symbols{i}, lines{i}, 1, 0);
        
        DT_data(index,:) = py';
        
        %Compute DTI from spline fit
        DTI(i,index) = 1 - (pmin.y)/(pmax.y);
        
        index = index + 1;
    end
    
    subplot(2, 2, ((i*2) - 1));
    contourf(px,time_slice,DT_data);
    colorbar;
    str = sprintf('%s: speed = %6.2f', FILE, unique_speed(i));
    title(str);
    
    DT_store{i} = DT_data;
    
    %get a resp-time slice at the best disparity
    [max_resp, max_disp_indx] = max(sum(DT_data));  %sum over time, take max over disparity
    PSTH_peak{i} = DT_data(:,max_disp_indx)-null_rate;
    PSTH_max(i) = max(PSTH_peak{i});
end


norm_max = max(PSTH_max); 
for i=1:length(unique_speed)	
    %normalize each slice to maximum of the pair (moving + static)
    PSTH_peak{i} = PSTH_peak{i}/norm_max;
end

subplot(2, 2, 2);
hold on;
plot(time_slice', DDI(1,:), 'r--');
plot(time_slice', PSTH_peak{1}, 'k--');
if (length(unique_speed) > 1)
    plot(time_slice', DDI(2,:), 'r-');
    plot(time_slice', PSTH_peak{2}, 'k-');
end
hold off;
ylim([0 1]);
%str = sprintf('%s: speed = %6.2f', FILE, unique_speed(i));
%title(str);

fileid = [BASE_PATH 'ProtocolSpecific\HDispTuning\PopTimeEvolution.mat'];
if exist(fileid, 'file') ~= 0
    load (fileid)
else 
    num_files = 0;
    stat_DTI = [];
    stat_DDI = [];
    stat_resp = [];
    mov_DTI = [];
    mov_DDI = [];
    mov_resp = [];
    sites = [];
    time = [];
end   

if (length(unique_speed) == 2)
    stat_DTI = [stat_DTI;  DTI(1,:) ];      
    stat_DDI = [stat_DDI;  DDI(1,:) ];      
    stat_resp = [stat_resp;  PSTH_peak{1}' ];
    mov_DTI = [mov_DTI;  DTI(2,:) ];      
    mov_DDI = [mov_DDI;  DDI(2,:) ];      
    mov_resp = [mov_resp;  PSTH_peak{2}' ];
end
if ( (length(unique_speed) == 1) & (unique_speed(1) == 0) )
    stat_DTI = [stat_DTI;  DTI(1,:) ];      
    stat_DDI = [stat_DDI;  DDI(1,:) ];      
    stat_resp = [stat_resp;  PSTH_peak{1}' ];
end
if ( (length(unique_speed) == 1) & (unique_speed(1) ~= 0) )
    mov_DTI = [mov_DTI;  DTI(1,:) ];      
    mov_DDI = [mov_DDI;  DDI(1,:) ];      
    mov_resp = [mov_resp;  PSTH_peak{1}' ];
end

time = [time; time_slice];

num_files = num_files + 1;
sites{num_files} = FILE;      
save (fileid, 'stat_DTI', 'stat_DDI', 'stat_resp', 'mov_DTI', 'mov_DDI', 'mov_resp', 'num_files', 'sites', 'time');


%yl = YLim;
%YLim([0 yl(2)]);	% set the lower limit of the Y axis to zero
%XLabel('Horizontal Disparity(deg)');
%YLabel('Response (spikes/sec)');

%now, print out some useful information in the upper subplot
%subplot(2, 1, 1);
%PrintGeneralData(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

return;