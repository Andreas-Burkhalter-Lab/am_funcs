%-----------------------------------------------------------------------------------------------------------------------
%-- HDispTuning_Playback.m -- Play back data from a horizontal disparity tuning run
%--	GCD, 4/22/01
%-----------------------------------------------------------------------------------------------------------------------
function HDispTuning_Playback(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);

TEMPO_Defs;		%needed for defines like IN_T1_WIN_CD
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01

        HDispTuning_Playback_Loop(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);%, gui_handles);
return;


%----------------------------------------------------------------------------------------------------------------------------------------------------
function HDispTuning_Playback_Loop(data, Protocol, Analysis, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE);%, gui_handles);

TEMPO_Defs;		%needed for defines like IN_T1_WIN_CD
ProtocolDefs;	%needed for all protocol specific functions - contains keywords - BJP 1/4/01

%get the column of values of disparities in the dots_params matrix
hdisp = data.dots_params(DOTS_HDISP,:,PATCH1);
unique_hdisp = munique(hdisp');

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
Pref_HDisp = data.one_time_params(PREFERRED_HDISP);

figure(2);
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [2 55 900 600], 'Name', 'Disparity Tuning');

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
select_trials = find( (hdisp >= -20) & (hdisp <= 20) );

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
    h_spikes = subplot('position', [.07 .77 .86 .2]);
    subplot(h_spikes);
    xtim = 1:length(spikes);
    h_raster = plot((xtim-data.htb_header{SPIKE_DB}.offset)*spike_bin_width, spikes, 'w-', 'LineWidth', 1.3);
    set(h_raster, 'EraseMode', 'none');
    hold on;
    plot([0 stim_duration], [1.2 1.2], 'k-', 'LineWidth', 4);
    hold off;
    xlim([-1 3]);
    ylim([0 1.5])
    xlabel('Time (seconds)', 'FontSize', 16);
    pause(0.5);

    %widen the width of each spike pulse, to make them more audible
    ii = find(spikes == 1);
    spikes(ii + 1) = 1;
    spikes(ii + 2) = 1;
    spikes(ii + 3) = 1;
        
    %convert to audio sampling rate, 8192 Hz
    spikes_interp = interp(spikes, 8); %this is an approximation
    sound(8*spikes_interp);

    h_eyes = subplot('position', [.07 .08 .4 .6]);
    subplot(h_eyes);

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
     
    hold on;
    h_cursor = plot(data.eye_data(LEYE_H,1,select_trials(i)), data.eye_data(LEYE_V,1,select_trials(i)), 'r+', 'MarkerSize', 10, 'LineWidth', 2);
    set(h_cursor, 'EraseMode', 'xor');

    h_dots = plot(dots_x_start, dots_y_start, 'wo', 'MarkerSize', 3, 'MarkerFaceColor', 'w');
    set(h_dots, 'EraseMode', 'xor');
    
    %plot a circle to show the RF location
    rectangle('Curvature', [1 1], 'Position',[(rf_x-rf_diam/2) (rf_y-rf_diam/2) rf_diam rf_diam]);
    %plot boxes to show the FP and targ windows
    for j=FP:FP
        rectangle('Position',[(targ_x(j)-targ_wid(j)/2) (targ_y(j)-targ_hgt(j)/2) targ_wid(j) targ_hgt(j)]);
    end
    %plot a marker for the fixation point
    plot(targ_x(FP), targ_y(FP), 'ks', 'MarkerSize', 4, 'MarkerFaceColor', 'k');
    
    %set the axis limits
    xlim([-10 10]);
    ylim([-10 10]);
    xlabel('Horizontal Position (deg)', 'FontSize', 16);
    ylabel('Vertical Position (deg)', 'FontSize', 16);

    count = 1;
    while (elap < 3.9)

        %compute an index into the spike array
        spk_index = floor(elap/spike_bin_width + 1);
        if ( mod(count,4) == 0 )
            set(h_raster, 'XData', (xtim(1:spk_index)-data.htb_header{SPIKE_DB}.offset)*spike_bin_width, 'YData', spikes(1:spk_index), 'Color', 'b' );
        end
        
        %compute an index into the eye position array
        index = floor(elap/eye_bin_width + 1);
        %plot the eye cursor
        if (data.eye_data(LEYE_H,index,select_trials(i)) ~= 0)
            set(h_cursor, 'XData', data.eye_data(LEYE_H,index,select_trials(i)), 'YData', data.eye_data(LEYE_V,index,select_trials(i)));
        end
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
            set(h_dots, 'XData', dots_x, 'YData', dots_y, 'MarkerEdgeColor', 'k', 'MarkerSize', 3, 'MarkerFaceColor', 'k' );
            
            buff = sprintf('Disparity = %6.3f deg', hdisp(select_trials(i)) );
            text(0, 6, buff, 'FontSize', 18, 'HorizontalAlignment', 'Center');
        end
        
        hold off;
        %force the figure to redraw
        drawnow;

        %get the current elapsed time
        elap = etime(clock, t0);
        count = count + 1;
    end
    cla;    
    
    %update the tuning curve
    h_curve = subplot('position', [.55 .08 .4 .6]);
    subplot(h_curve);    
    [xx, index] = sort(hdisp(select_trials(1:i)));
    yy = spike_rates(select_trials(1:i));
    yy = yy(index);
    %[xx' yy']
    PlotTuningCurve(xx, yy, 'ko', 'k-', 0, 1);
    xlim([-2 2]);
    ylim([0 max(spike_rates)]);
    xlabel('Horizontal Disparity (deg)', 'FontSize', 16);
    ylabel('Response (spikes/sec)', 'FontSize', 16);
    
    subplot(h_spikes);
    cla;
    pause(1);
end

return;