%-----------------------------------------------------------------------------------------------------------------------
%-- SpeedTuning_Playback.m -- Play back data from a speed tuning run
%--	GCD, 4/22/01
%-----------------------------------------------------------------------------------------------------------------------
function SpeedTuning_Playback(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;		%needed for defines like IN_T1_WIN_CD
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01


%global status_flag;

%gui = Playback_GUI;
%gui_handles = guihandles(gui);

%set(gui_handles.pushbutton1, 'Value', 1);

%while(1)
%    pause(1);
    %temp = get(gui_handles.pushbutton1, 'Value')
    %if (temp == 1)
    %status_flag
    %if (status_flag == 1)
        SpeedTuning_Playback_Loop(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);%, gui_handles);
        %set(gui_handles.pushbutton1, 'Enable', 'on')
        %set(gui_handles.pushbutton2, 'Enable', 'off')
        %end
        %end
return;


%----------------------------------------------------------------------------------------------------------------------------------------------------
function SpeedTuning_Playback_Loop(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);%, gui_handles);

TEMPO_Defs;		%needed for defines like IN_T1_WIN_CD
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01

%get the column of values of directions in the dots_params matrix
direction = data.dots_params(DOTS_DIREC,:,PATCH1);
unique_direc = munique(direction');

%get the column of values of speeds in the dots_params matrix
speed = data.dots_params(DOTS_SPEED,:,PATCH1);
unique_speed = munique(speed');

%now, get the firing rates for all the trials 
spike_rates = data.spike_rates(SpikeChan, :);

Pref_Direc = data.one_time_params(PREFERRED_DIRECTION);
Pref_Speed = data.one_time_params(PREFERRED_SPEED);

figure(2);
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [485 440 540 298], 'Name', 'Spike Train','MenuBar', 'None');
figure(3);
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [485 35 540 380], 'Name', 'Speed Tuning Curve','MenuBar', 'None');
figure(4);
%set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [2 35 475 475], 'Name', 'Behavior Window','MenuBar', 'None');
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [2 35 475 475], 'Name', 'Behavior Window');

rf_x = data.one_time_params(RF_XCTR);
rf_y = data.one_time_params(RF_YCTR);
rf_diam = data.one_time_params(RF_DIAMETER);

%first, we'll need the bin_width (in sec) of our spike raster + events + eye data
h = data.htb_header{SPIKE_DB};
spike_bin_width = (h.skip + 1) / (h.speed_units / h.speed);
h = data.htb_header{EVENT_DB};
event_bin_width = (h.skip + 1) / (h.speed_units / h.speed);
h = data.htb_header{EYE_DB};
eye_bin_width = (h.skip + 1) / (h.speed_units / h.speed);

%pull out trials for the selected correlation level
select_trials = find(speed >= 0);

%loop through all the selected trials
for i = 1:length(select_trials)

    %find event bins corresponding to STIM_START and STIM_STOP
    start_eventbin = find(data.event_data(1,:,select_trials(i)) == VSTIM_ON_CD);
    stop_eventbin = find(data.event_data(1,:,select_trials(i)) == VSTIM_OFF_CD);
    stim_duration = (stop_eventbin - start_eventbin)*event_bin_width;
    
    spikes = [];
    spikes = data.spike_data(SpikeChan,:,select_trials(i));
    bad_ones = (spikes > 1); %shouldn't happen, but just in case
    spikes(logical(bad_ones)) = 1;
    
    %plot the spike train
    figure(2);
    xtim = 1:length(spikes);
    plot((xtim-data.htb_header{SPIKE_DB}.offset)*spike_bin_width, spikes, 'b-', 'LineWidth', 1.3);
    hold on;
    plot([0 stim_duration], [1.2 1.2], 'k-', 'LineWidth', 4);
    hold off;
    ylim([0 1.5])
    xlabel('Time (seconds)');
    pause(0.5);

    %widen the width of each spike pulse, to make them more audible
    ii = find(spikes == 1);
    spikes(ii + 1) = 1;
    spikes(ii + 2) = 1;
    spikes(ii + 3) = 1;
        
    %convert to audio sampling rate, 8192 Hz
    spikes_interp = interp(spikes, 8); %this is an approximation
    %max(spikes_interp)
    sound(8*spikes_interp);
    figure(4);
    %set(gcf, 'Color', 'White');
    count = 0;
    t0 = clock;
    elap = etime(clock, t0);
    targ_x =   [data.targ_params(TARG_XCTR,select_trials(i),FP) data.targ_params(TARG_XCTR,select_trials(i),T1) data.targ_params(TARG_XCTR,select_trials(i),T2)];
    targ_y =   [data.targ_params(TARG_YCTR,select_trials(i),FP) data.targ_params(TARG_YCTR,select_trials(i),T1) data.targ_params(TARG_YCTR,select_trials(i),T2)];
    targ_wid = [data.targ_params(TARG_WIN_FWIDTH,select_trials(i),FP) data.targ_params(TARG_WIN_FWIDTH,select_trials(i),T1) data.targ_params(TARG_WIN_FWIDTH,select_trials(i),T2)];
    targ_hgt = [data.targ_params(TARG_WIN_FHEIGHT,select_trials(i),FP) data.targ_params(TARG_WIN_FHEIGHT,select_trials(i),T1) data.targ_params(TARG_WIN_FHEIGHT,select_trials(i),T2)];
    
    dots_x_start = []; dots_y_start = []; dots_x = []; dots_y = [];
    % compute some random dot locations
    num_dots = 50;
    dots_x_start = rand(num_dots, 1)*rf_diam + rf_x - rf_diam/2;
    dots_y_start = rand(num_dots, 1)*rf_diam + rf_y - rf_diam/2;
     
    while (elap < 3.9)
        %compute an index into the eye position array
        index = floor(elap/eye_bin_width + 1);
        %plot the eye cursor
        if (data.eye_data(LEYE_H,index,select_trials(i)) ~= 0)
            plot(data.eye_data(LEYE_H,index,select_trials(i)), data.eye_data(LEYE_V,index,select_trials(i)), 'r+', 'MarkerSize', 10, 'LineWidth', 2);
        end
        hold on;
        %plot boxes to show the FP and targ windows
        for j=FP:FP
            rectangle('Position',[(targ_x(j)-targ_wid(j)/2) (targ_y(j)-targ_hgt(j)/2) targ_wid(j) targ_hgt(j)]);
        end
        %plot a circle to show the RF location
        rectangle('Curvature', [1 1], 'Position',[(rf_x-rf_diam/2) (rf_y-rf_diam/2) rf_diam rf_diam]);
        %plot a marker for the fixation point
        plot(targ_x(FP), targ_y(FP), 'ks', 'MarkerSize', 4, 'MarkerFaceColor', 'k');
        %during the visual stimulus period, plot text in RF to show stimulus polarity
        if ( (elap > start_eventbin*event_bin_width) & (elap < stop_eventbin*event_bin_width) )  %the visual stimulus period
            dots_x = dots_x_start + speed(select_trials(i))*(elap-start_eventbin*event_bin_width)*cos(direction(select_trials(i))*pi/180);
            while ( sum((dots_x>(rf_x+rf_diam/2))) > 0)
                dots_x(dots_x>(rf_x+rf_diam/2)) = dots_x(dots_x>(rf_x+rf_diam/2)) - rf_diam;
            end
            while ( sum((dots_x<(rf_x-rf_diam/2))) > 0)
                dots_x(dots_x<(rf_x-rf_diam/2)) = dots_x(dots_x<(rf_x-rf_diam/2)) + rf_diam;
            end

            dots_y = dots_y_start + speed(select_trials(i))*(elap-start_eventbin*event_bin_width)*sin(direction(select_trials(i))*pi/180);
            while ( sum((dots_y>(rf_y+rf_diam/2))) > 0)
                dots_y(dots_y>(rf_y+rf_diam/2)) = dots_y(dots_y>(rf_y+rf_diam/2)) - rf_diam;
            end
            while ( sum((dots_y<(rf_y-rf_diam/2))) > 0)
                dots_y(dots_y<(rf_y-rf_diam/2)) = dots_y(dots_y<(rf_y-rf_diam/2)) + rf_diam;
            end

            radius = sqrt((dots_x - rf_x).^2 + (dots_y - rf_y).^2);
            dots_x(radius > rf_diam/2) = NaN;
            dots_y(radius > rf_diam/2) = NaN;
            plot(dots_x, dots_y, 'ko', 'MarkerSize', 3, 'MarkerFaceColor', 'k');
            
            buff = sprintf('Speed = %5.1f deg per sec', speed(select_trials(i)) );
            text(0, 9, buff, 'FontSize', 18, 'HorizontalAlignment', 'Center');
        end
        
        %set the axis limits
        xlim([-20 20]);
        ylim([-20 20]);
        hold off;
        %force the figure to redraw
        drawnow;
        %get the current elapsed time
        elap = etime(clock, t0);
    end
    
    
    %update the response histogram
    figure(3);
    clf;
    [xx, index] = sort(speed(select_trials(1:i)));
    yy = spike_rates(select_trials(1:i));
    yy = yy(index);
    %[xx' yy']
    PlotTuningCurve(xx, yy, 'ko', 'k-', 0, 1);
    xlim([-5 35]);
    ylim([0 max(spike_rates)]);
    xlabel('Speed (deg per sec)');
    ylabel('Response (spikes/sec)');
    
    figure(2);
    clf;
    pause(1);

    %if (get(gui_handles.pushbutton1, 'Value') == 0)    
    %if (status_flag == 0)    
    %    break;
    %end
end


return;