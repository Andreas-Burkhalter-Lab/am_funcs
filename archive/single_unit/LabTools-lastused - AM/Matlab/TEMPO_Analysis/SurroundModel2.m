function SurroundModel2(data, SpikeChan, StartCode, StopCode, BegTrial, EndTrial, StartOffset, StopOffset, PATH, FILE, HDisparity, SurrMap, UseSyncPulses);

TEMPO_Defs;

%following function checks offset times and outputs starting and ending analysis (in spike bins)
num_trials = size(data.event_data, 3);
[StartOffsetBin StopOffsetBin] = CheckTimeOffset(data, num_trials, StartCode, StopCode, StartOffset, StopOffset, UseSyncPulses);

if (~isempty(data.spike_data))
    %compute the firing rate over all trials during the period between StartCode and StopCode
    data.spike_rates = ComputeSpikeRates(data, num_trials, StartCode, StopCode, StartOffsetBin, StopOffsetBin);
end

symbols = {'bo' 'ro' 'go' 'ko' 'b*' 'r*' 'g*' 'k*'};
lines = {'b-' 'r-' 'g-' 'k-' 'b--' 'r--' 'g--' 'k--'};

%get position of peak response within the size tuning curve
%use this as the size of the receptive field (will be updated later for more sophisticated
%size analyses.

%get the column of values of horiz. disparity in the dots_params matrix
hor_disp = data.dots_params(DOTS_HDISP,BegTrial:EndTrial,PATCH1);

%get the column of values of horiz. disparity magnitude in the dots_params matrix
mag_disp = data.dots_params(DOTS_HGRAD_MAG,BegTrial:EndTrial,PATCH1);

%get indices of any NULL conditions (for measuring spontaneous activity)
null_trials = logical((mag_disp == data.one_time_params(NULL_VALUE)) );

unique_mag_disp = munique(mag_disp(~null_trials)');	

%get the column of values of horiz. disparity angle of orientation in the dots_params matrix
disp_ang = data.dots_params(DOTS_HGRAD_ANGLE,BegTrial:EndTrial,PATCH1);
unique_disp_ang = munique(disp_ang(~null_trials)');

%get the column of mean disparity values
mean_disp = data.dots_params(DOTS_HDISP,BegTrial:EndTrial,PATCH1);
unique_mean_disp = munique(mean_disp(~null_trials)');

%get the column of values of directions in the dots_params matrix
ap_size = data.dots_params(DOTS_AP_XSIZ,:,PATCH1);
unique_size = munique(ap_size(~null_trials)');

%now, get the firing rates for all the trials 
spike_rates = data.spike_rates(SpikeChan, :);

%get indices of any NULL conditions (for measuring spontaneous activity
null_trials = logical( (ap_size == data.one_time_params(NULL_VALUE)) );

%get indices of monoc. and uncorrelated controls
control_trials = logical( (mag_disp == LEYE_CONTROL) | (mag_disp == REYE_CONTROL) | (mag_disp == UNCORR_CONTROL) );

%now, remove trials from direction and spike_rates that do not fall between BegTrial and EndTrial
trials = 1:length(ap_size);		% a vector of trial indices
select_trials = ( (trials >= BegTrial) & (trials <= EndTrial) );

%this will be the diameter of the grid
CRF_Size = SurrMap.ctrsize;

Total_RF_Size = SurrMap.ctrsize + (SurrMap.surrsize * 2);
Surround_Size = SurrMap.surrsize;
%number of surround patches
num_surr = length(SurrMap.angles);

%retrieve the tuning curves for the CRF and the surrounds
CRF.px = HDisparity.px;
CRF.py = HDisparity.py;
CRF.plot_x = HDisparity.plot_x;
CRF.plot_y = HDisparity.plot_y;

%pad with uncorrelated values
%CRF.uncorr = HDisparity.uncorr_resp;
%CRF.px = [-CRF_Size; CRF.px; CRF_Size];
%CRF.py = [CRF.uncorr; CRF.py; CRF.uncorr];
%CRF.plot_x = [-CRF_Size CRF.plot_x CRF_Size];
%CRF.plot_y = [CRF.uncorr CRF.plot_y CRF.uncorr];

for i = 1:num_surr
    temp_px = SurrMap.px(:, i);
    temp_py = SurrMap.py(:, i);
    temp_plot_x = SurrMap.plot_x(i, :);
    temp_plot_y = SurrMap.plot_y(i, :);

    %SurrPatch.uncorr(:,i) = HDisparity.uncorr_resp;
    %SurrPatch.px(:,i) = [-Surround_Size; temp_px; Surround_Size];
    %SurrPatch.py(:,i) = [SurrPatch.uncorr(:,i); temp_py; SurrPatch.uncorr(:,i)];    
    %SurrPatch.plot_x(i, :) = [-Surround_Size temp_plot_x Surround_Size];
    %SurrPatch.plot_y(i, :) = [SurrPatch.uncorr(:,i) temp_plot_y SurrPatch.uncorr(:,i)];
    SurrPatch.px(:,i) = SurrMap.px(:, i);           %[-Surround_Size; temp_px; Surround_Size];
    SurrPatch.py(:,i) = SurrMap.py(:, i);           %[SurrPatch.uncorr(:,i); temp_py; SurrPatch.uncorr(:,i)];    
    SurrPatch.plot_x(i, :) = SurrMap.plot_x(i, :);  %[-Surround_Size temp_plot_x Surround_Size];
    SurrPatch.plot_y(i, :) = SurrMap.plot_y(i, :);  %[SurrPatch.uncorr(:,i) temp_plot_y SurrPatch.uncorr(:,i)];
end

%figure that displays tuning curves in appropriate locations
surr_tuning_fig=figure;

crf_axis = axes('position', [.35 .32 .3 .3]);
surr_axis(1) = axes('position', [.68 .32 .3 .3]);
surr_axis(2) = axes('position', [.52 .64 .3 .3]);
surr_axis(3) = axes('position', [.18 .64 .3 .3]);
surr_axis(4) = axes('position', [.02 .32 .3 .3]);
surr_axis(5) = axes('position', [.18 0 .3 .3]);
surr_axis(6) = axes('position', [.52 0 .3 .3]);


%get ctr only response from surround tuning data
CTR_only = SurrMap.ctronly;

dogabor = 1;
dosurroundgabor = 1;

if dogabor == 1
    means = [CRF.px CRF.py];
    raw = [CRF.plot_x' CRF.plot_y'];
    fixed_param_flags = zeros(6,1); %by default, all 6 parameters will vary
    fixed_param_values = zeros(6,1); 
    [pars,freq] = gaborfit_for_SurroundModel(means,raw, fixed_param_flags, fixed_param_values, HDisparity.uncorr_resp);
    CRF.pars = pars;
    CRF.x_interp = (CRF.px(1)): .01 : (CRF.px(length(CRF.px)));
    % Note: this func MUST match that in gaborfit.m
    CRF.y_interp =  pars(1) + pars(2)*exp(-0.5*((CRF.x_interp - pars(3))/ pars(4)).^2).*cos(2*pi*pars(5)*(CRF.x_interp - pars(3))+pars(6) );
    CRF.y_interp(CRF.y_interp < 0) = 0;
    
    %scale response to account for small size of center patch used surround
    %model protocol
    Response_of_CRF_to_disp_in_ctr_patch = pars(1) + pars(2)*exp(-0.5*((SurrMap.ctrdisp - pars(3))/ pars(4)).^2).*cos(2*pi*pars(5)*(SurrMap.ctrdisp - pars(3))+pars(6) );
    scale_responses = Response_of_CRF_to_disp_in_ctr_patch/CTR_only;
    CRF.y_interp = CRF.y_interp/scale_responses;
    
else
    CRF.x_interp = (CRF.px(1)): .01 : (CRF.px(length(CRF.px)));
    CRF.y_interp = interp1(CRF.px, CRF.py, CRF.x_interp, 'spline');  %spline fit
    CRF_max = max(CRF.y_interp);
    CRF_diff = CRF_max - CTR_only;
    CRF.y_interp = CRF.y_interp - CRF_diff;
end

%tuning = figure;
%hold on
%plot(CRF.x_interp, CRF.y_interp, 'k--');
%plot crf tuning curve in center
figure(surr_tuning_fig)
axes(crf_axis)
plot(CRF.x_interp, CRF.y_interp, 'k--');

for i=1:num_surr
    if dosurroundgabor == 1
        means = [SurrPatch.px(:, i) SurrPatch.py(:, i)];
        raw = [SurrPatch.plot_x(i, :)' SurrPatch.plot_y(i, :)'];
        fixed_param_flags = zeros(6,1); %by default, all 6 parameters will vary
        fixed_param_values = zeros(6,1); 
        [pars,freq] = gaborfit_for_SurroundModel(means,raw, fixed_param_flags, fixed_param_values, HDisparity.uncorr_resp);
        SurrPatch.pars(:,i) = pars;
        SurrPatch.x_interp(i,:) = (SurrPatch.px(1,i)): .01 : (SurrPatch.px(length(SurrPatch.px(:,i)), i));
        SurrPatch.y_interp(i,:) =  pars(1) + pars(2)*exp(-0.5*((SurrPatch.x_interp(i,:) - pars(3))/ pars(4)).^2).*cos(2*pi*pars(5)*(SurrPatch.x_interp(i,:) - pars(3))+pars(6) );
        SurrPatch.y_interp(SurrPatch.y_interp(i,:) < 0) = 0;
    else
        SurrPatch.x_interp(i,:) = (SurrPatch.px(1,i)): .01 : (SurrPatch.px(length(SurrPatch.px(:,i)), i));
        SurrPatch.y_interp(i,:) = interp1(SurrPatch.px(:, i), SurrPatch.py(:, i), SurrPatch.x_interp(i,:), 'spline') - CTR_only;  %spline fit
    end
    %figure(tuning)
    %hold on
    %plot(SurrPatch.x_interp(i,:),SurrPatch.y_interp(i,:), lines{i});
end

for i=1:num_surr
    pars = SurrPatch.pars(:,i);
    x_interp_temp = (SurrPatch.px(1,i))-2: .01 : (SurrPatch.px(length(SurrPatch.px(:,i)), i))+2;
    y_interp_temp = pars(1) + pars(2)*exp(-0.5*((x_interp_temp - pars(3))/pars(4)).^2).*cos(2*pi*pars(5)*(x_interp_temp-pars(3))+pars(6));
    y_interp_gauss = pars(1) + pars(2)*exp(-0.5*((x_interp_temp - pars(3))/pars(4)).^2);
    figure(surr_tuning_fig);
    hold on
    axes(surr_axis(i));
    PlotTuningCurve(SurrPatch.plot_x(i, :), SurrPatch.plot_y(i, :), symbols{i}, '', 0, 1);
    hold on
    plot([SurrPatch.px(1,i)-2 (SurrPatch.px(length(SurrPatch.px(:,i)), i))+2], [HDisparity.uncorr_resp HDisparity.uncorr_resp], 'b--');
    hold on
    plot(x_interp_temp, y_interp_temp, lines{i});
    plot(x_interp_temp, y_interp_gauss, 'k--');
end

print_pars = 1;
if print_pars == 1
    filesize = size(FILE,2) - 1;
    while FILE(filesize) ~='.'
        filesize = filesize - 1;
    end
    PATHOUT = 'Z:\Users\jerry\GradAnalysis\';
    outfile = [PATHOUT 'Surround_parameters_for_' FILE(1:filesize) 'dat'];
    fid = fopen(outfile, 'a');
    for i = 1:num_surr
        pars = SurrPatch.pars(:,i);
        fprintf(fid, '%3.2f\t%3.2f\t%3.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\n', SurrMap.angles(i), pars(1), pars(2), pars(3), pars(4), pars(5), pars(6)); 
    end
    fclose(fid);
end

%figure(tuning)
%angle1 = sprintf('%d', SurrMap.angles(1));
%angle2 = sprintf('%d', SurrMap.angles(2));
%angle3 = sprintf('%d', SurrMap.angles(3));
%angle4 = sprintf('%d', SurrMap.angles(4));
%angle5 = sprintf('%d', SurrMap.angles(5));
%angle6 = sprintf('%d', SurrMap.angles(6));
%legend('CRF', angle1, angle2, angle3, angle4, angle5, angle6, 0);

surr_max = max(SurrPatch.px(:,1));
surr_min = min(SurrPatch.px(:,1));

%draw a little schematic of what the surround and center look like
CRF_xdraw = 0-(CRF_Size/2);
CRF_ydraw = 0-(CRF_Size/2);
%circles = figure;
%hold on
%rectangle('Position', [-unique_size(1)/2 -unique_size(1)/2 unique_size(1) unique_size(1)], 'Curvature', [1 1], 'FaceColor', 'b');
%rectangle('Position', [CRF_xdraw CRF_ydraw CRF_Size CRF_Size], 'Curvature', [1 1], 'FaceColor','r');
%axis equal

Distance_between_Surr_and_Ctr_ctrs = (Surround_Size + CRF_Size)/2;

surr_xdraw = [];
surr_ydraw = [];
for i=1:num_surr
    ang = SurrMap.angles(i);
    ang_rad = ang * (pi/180);
    surr_xdraw(i) = Distance_between_Surr_and_Ctr_ctrs * cos(ang_rad) + 0 - Surround_Size/2;
    surr_ydraw(i) = Distance_between_Surr_and_Ctr_ctrs * sin(ang_rad) + 0 - Surround_Size/2;
    %  rectangle('Position', [surr_xdraw(i) surr_ydraw(i) Surround_Size Surround_Size], 'Curvature', [1 1]);
end
%resolution of grids in degrees
res = .1;

x_offset = abs(min(surr_xdraw)) + 1;
y_offset = abs(min(surr_ydraw)) + 1;

x_ind = fix((surr_xdraw + x_offset) * 1/res);
y_ind = fix((surr_ydraw + y_offset) * 1/res);

x_max = max(x_ind) + Surround_Size*(1/res);
y_max = max(y_ind) + Surround_Size*(1/res);

CRF_xind = fix((CRF_xdraw + x_offset) * 1/res);
CRF_yind = fix((CRF_ydraw + y_offset) * 1/res);

zero_xind = fix((0+x_offset) * 1/res);
zero_yind = fix((0+y_offset) * 1/res);


%big_circles = figure;
%hold on
%rectangle('Position', [-unique_size(1)/2 -unique_size(1)/2 unique_size(1) unique_size(1)], 'Curvature', [1 1], 'FaceColor', 'b');
%rectangle('Position', [CRF_xind CRF_yind CRF_Size*(1/res) CRF_Size*(1/res)], 'Curvature', [1 1], 'FaceColor','r');
%axis equal
%for i=1:num_surr
%   rectangle('Position', [x_ind(i) y_ind(i) Surround_Size*(1/res) Surround_Size*(1/res)], 'Curvature', [1 1]);
%end
%axis([0 x_max 0 y_max]);

%grid_size = max([x_max y_max]);
x_list = 1:x_max;
y_list = 1:y_max;
x_coord = (x_list-zero_xind) * res;
y_coord = (y_list-zero_yind) * res;

distance_grid = zeros(x_max, y_max);

for i = 1:length(x_coord)
    for j = 1:length(y_coord)
        distance_grid(i, j) = sqrt((x_coord(i)-0)^2 + (y_coord(j)-0)^2);
    end
end

grad_fig = figure;

disparity_and_angles = zeros(8, length(unique_disp_ang));

pref_tilt = zeros(1, length(unique_mean_disp));
pref_tilt_meas = zeros(1, length(unique_mean_disp));
p_val = zeros(1, length(unique_mean_disp));

%figure to contains the value of the mean disp and model value for each
%position
tilt_text_fig = figure;
set(gcf,'PaperPosition', [.2 .2 8 10.7], 'Position', [500 50 800 773]);
subplot(3,2,1);
%surr_axis2(1) = axes('position', [.68 .32 .5 .3]);
%surr_axis2(2) = axes('position', [.68 .64 .5 .3]);
%surr_axis2(3) = axes('position', [.02 .64 .5 .3]);
%surr_axis2(4) = axes('position', [.02 .32 .5 .3]);
%surr_axis2(5) = axes('position', [.02 0 .5 .3]);
%surr_axis2(6) = axes('position', [.68 0 .5 .3]);

angle_out = '';
for magdisp =1:length(unique_mag_disp) %for each magnitude of disparity
    resp_mag_ang = zeros(length(unique_mag_disp), length(unique_disp_ang));
    all_data = [];
    angle_data = [];
    grad_mean = zeros(length(unique_disp_ang), 2);
    data_string = [];
    m_disp_max = zeros(length(unique_mean_disp), 1);
    m_disp_min = zeros(length(unique_mean_disp), 1);
    
    %variables for TDI calculation
    start = zeros(length(unique_mean_disp), 1);
    stop = zeros(length(unique_mean_disp), 1);
    start_save = zeros(length(unique_mean_disp), 1);
    stop_save = zeros(length(unique_mean_disp), 1);
    TDIdata = [];
    Save_Responses = [];
    ax_x = -.5;
    for mdisp=1:length(unique_mean_disp)  %for each mean disparity...
        number_disp_outside_range_per_meandisp = 0;
        angle_and_response_data = [];
        ax_x = ax_x + .4;
        ax_y = -.05;
        ax_y = ax_y + .05;
        for i=1:num_surr
            figure(tilt_text_fig)
            if i == 1
                subplot(3, 2, 4)
            elseif i == 2
                subplot(3,2,2);
            elseif i == 3
                subplot(3,2,1);
            elseif i == 4
                subplot(3,2,3);
            elseif i == 5
                subplot(3,2,5);
            elseif i == 6
                subplot(3,2,6);
            end
            hold on
            axis ij
            axis off
            string_out = sprintf('%3.2f', unique_mean_disp(mdisp));
            text(ax_x, ax_y, string_out, 'FontSize', 8);                
        end
        for dispang =1:length(unique_disp_ang)

            %create a matrix where each cell contains the disparity of that location
            disparity_grid = zeros(x_max, y_max) + unique_mean_disp(mdisp);
            
            angle = unique_disp_ang(dispang) * (pi/180);
            
            for l = 1:x_max %for each x pt on the grid
                %get the x coordinate of the grid
                x_val = x_coord(l);
                for m = 1:y_max %for each y point on the grid
                    %get the y coordinate of the grid
                    y_val = y_coord(m);
                    if(angle==0)
                        disp_val = y_val * unique_mag_disp(magdisp);
                        %is alpha != 0 then there will be a tilt to the gradient
                    else
                        len = distance_grid(l,m);
                        if len ~= 0
                            theta = asin(y_val/len);
                            if ~isreal(theta)
                                theta = real(theta);
                            end
                            %//sin only goes from -90 to 90 degrees.  need to fix for the other two quadrants.
                            %//need to do quadrant corrections
                            %//cases (-, +), (-, -)
                            if ((x_val <= 0) & (y_val >= 0))
                                %//Quadrant III
                                %//theta is now greater than 90 degrees
                                temptheta = (pi/2)-theta;
                                theta = (pi/2) + temptheta;
                            elseif ((x_val <= 0) & (y_val <= 0))
                                %//Quadrant IV
                                %//theta is now less than -90 degrees
                                temptheta = (-pi/2)-theta;
                                theta = (-pi/2) + temptheta;
                            end %end quad check if
                            
                            %//calculate the distance of the point from the gradient tilt axis
                            y_tilt = len * sin(theta-angle);
                            
                            %//calculate the disparity difference
                            disp_val = y_tilt * unique_mag_disp(magdisp);
                        else
                            disp_val = 0;
                        end %end len if
                    end %end angle if
                    disparity_grid(l,m) = disparity_grid(l,m) + disp_val;
                end %end m
            end %end l
            %now that we have the disparity gradient we can calculate the mean disparity
            %in each of the surround patches
            mean_disp_in_surr = [];
            resp = [];
            response = 0;
            
            %switches for choosing which interpolation scheme to use
            dolinear = 1;
            docubic = 0;
            dospline = 0;
            
            number_disp_outside_range_per_angle = 0;
            
            %calculate mean disparity for each surround patch
            %then interpolate using the surround tuning curves to calculate the expected response of the surround
            %to that disparity
            %num_surr = 6;
            ax_y = ax_y + .1;
            for i=1:num_surr

                disp_in_surr = disparity_grid(x_ind(i):x_ind(i)+Surround_Size*(1/res), y_ind(i):y_ind(i)+Surround_Size*(1/res));
                mean_disp_in_surr(i) = mean(disp_in_surr(:));
                if mean_disp_in_surr(i)>surr_max | mean_disp_in_surr(i)<surr_min
                    number_disp_outside_range_per_angle = number_disp_outside_range_per_angle + 1;
                end
                
                %use model fit to calculate response from disparity
                pars = SurrPatch.pars(:,i);                
                resp(i) = pars(1) + pars(2)*exp(-0.5*((mean_disp_in_surr(i) - pars(3))/ pars(4)).^2).*cos(2*pi*pars(5)*(mean_disp_in_surr(i) - pars(3))+pars(6) );
                
                %print out resp and mean disp in the appropriate axis
                %location
                figure(tilt_text_fig);
                if i == 1
                    subplot(3, 2, 4)
                elseif i == 2
                    subplot(3,2,2);
                elseif i == 3
                    subplot(3,2,1);
                elseif i == 4
                    subplot(3,2,3);
                elseif i == 5
                    subplot(3,2,5);
                elseif i == 6
                    subplot(3,2,6);
                end
                hold on
                axis off
                axis ij
                string_print = sprintf('%3.2f, %1.3f, %3.2f', unique_disp_ang(dispang), mean_disp_in_surr(i), resp(i));
                text(ax_x, ax_y, string_print, 'FontSize', 8);
    
                if dolinear == 1
                    test(i) = interp1(SurrPatch.x_interp(i,:), SurrPatch.y_interp(i,:), mean_disp_in_surr(i), 'linear', 'extrap');
                elseif docubic == 1
                    test(i) = interp1(SurrPatch.x_interp(i,:), SurrPatch.y_interp(i,:), mean_disp_in_surr(i), 'cubic', 'extrap');
                elseif dospline == 1
                    test(i) = interp1(SurrPatch.x_interp(i,:), SurrPatch.y_interp(i,:), mean_disp_in_surr(i), 'spline', 'extrap');
                end
            end
            
            %calculate the mean disparity in the CRF
            %then using the disparity tuning curve of the center
            %inpterpolate to find the expected response
            mean_disp_and_res = [mean_disp_in_surr' resp'];
            disp_in_CRF = disparity_grid(CRF_xind:CRF_xind+CRF_Size*(1/res) ,  CRF_yind:CRF_yind + CRF_Size*(1/res));
            mean_disp_in_CRF = mean(disp_in_CRF(:));
            
            %use model fit to calculate response from disparity
            pars = CRF.pars;                
            CRF_response = pars(1) + pars(2)*exp(-0.5*((mean_disp_in_CRF - pars(3))/ pars(4)).^2).*cos(2*pi*pars(5)*(mean_disp_in_CRF - pars(3))+pars(6) );
            CRF_response = CRF_response/scale_responses;
            if dolinear == 1
                CRF_test = interp1(CRF.x_interp, CRF.y_interp, mean_disp_in_CRF, 'linear', 'extrap');
            elseif docubic == 1
                CRF_test = interp1(CRF.x_interp, CRF.y_interp, mean_disp_in_CRF, 'cubic', 'extrap');
            elseif dospline == 1
                CRF_test = interp1(CRF.x_interp, CRF.y_interp, mean_disp_in_CRF, 'spline', 'extrap');
            end
            
            disparities_and_angle(1,dispang) = unique_disp_ang(dispang);
            disparities_and_angle(2,dispang) = mean_disp_in_CRF;
            disparities_and_angle(3:8,dispang) = mean_disp_in_surr';
            
            %the "-6*CTR_only" part removes the response of the CRF from the responses of the test patches 
            %(since both were presented simultaneously            
            avg_response = CRF_response + sum(resp) - 6*CTR_only; 
            angle_and_response_data(dispang,:) = [unique_disp_ang(dispang) avg_response];
            
            number_disp_outside_range_per_meandisp = number_disp_outside_range_per_meandisp + number_disp_outside_range_per_angle;
        end %end for each angle
        
        figure(grad_fig)
        hold on
        
        %plot the real data with the model data
        sizem = length(unique_size);
        disp_select = logical((ap_size == unique_size(sizem)) & (mag_disp == unique_mag_disp(magdisp)) & (mean_disp == unique_mean_disp(mdisp)) );
        
        plot_x_measured = disp_ang(disp_select & ~null_trials & ~control_trials & select_trials);
        plot_y_measured = spike_rates(disp_select & ~null_trials & ~control_trials & select_trials);
        %NOTE: inputs to PlotTuningCurve must be column vectors, not row vectors, because of use of munique()
        %this is measured data plot
        [px_m, py_m, perr_m, spk_max_m, spk_min_m] = PlotTuningCurve(plot_x_measured', plot_y_measured', symbols{mdisp}, lines{mdisp}, 1, 0);
        
        %calculating max for measured--------------------------------------------------------
        %calculate max from curve fit
        px_measured_save = px_m;
        py_measured_save = py_m;
        px_m = (px_m * pi)/180;
        plot_x = (plot_x_measured * pi)/180;
        means = [px_m py_m];
        raw = [plot_x plot_y_measured];
        
        %fit with a distorted sin wave
        pars = sin_exp_fit_model(means,raw);
        
        x_interp = (px_m(1)): .01 : (px_m(length(px_m)));
        y_sin = sin_exp_func(x_interp, pars);
        y_err = sin_exp_err_model(pars);
        x_rad = (x_interp * 180)/pi;
        
        
        [max_y, ind] = max(y_sin);
        [min_y, ind_min] = min(y_sin);
        
        angle_max = x_rad(ind);
        angle_min = x_rad(ind_min);
        
        pref_tilt_meas(1, mdisp) = angle_max;
        
        %store p-values of each curve
        temp_x =plot_x_measured;
        p_val(1,mdisp) = calc_mdisp_anovap(disp_select, temp_x, plot_y_measured, unique_disp_ang);
        %end calculating max for measured--------------------------------------------------------
        
        %store pertinent information about measured data to measure TDI later
        start(mdisp) = length(TDIdata)+1;
        stop(mdisp) = length(plot_x_measured)+start(mdisp)-1;
        TDIdata(start(mdisp):stop(mdisp), 1) = plot_x_measured';
        TDIdata(start(mdisp):stop(mdisp), 2) = plot_y_measured';
        
        %hold on
        blanklines = '';
        [px, py, perr, spk_max, spk_min] = PlotTuningCurve(angle_and_response_data(:, 1), angle_and_response_data(:, 2), symbols{mdisp}, blanklines, 0, 1);
        legend_on = 1;
        
        px_model_save = px;
        py_model_save = py;
        
        %store pertinent information about responses (measured and model) for
        %correlation
        start_save(mdisp) = length(Save_Responses)+1;
        stop_save(mdisp) = length(px_measured_save)+start_save(mdisp)-1;
        Save_Responses(start_save(mdisp):stop_save(mdisp), 1) = px_measured_save;
        Save_Responses(start_save(mdisp):stop_save(mdisp), 2) = py_measured_save;
        Save_Responses(start_save(mdisp):stop_save(mdisp), 3) = px_model_save;
        Save_Responses(start_save(mdisp):stop_save(mdisp), 4) = py_model_save;
        
        %calculating max for model--------------------------------------------------------
        %calculate max from curve fit
        
        px = (px * pi)/180;
        plot_x = (angle_and_response_data(:, 1) * pi)/180;
        means = [px py];
        raw = [plot_x angle_and_response_data(:, 2)];
        
        %fit with a distorted sin wave
        pars = sin_exp_fit_model(means,raw);
        
        x_interp = (px(1)): .01 : (px(length(px)));
        y_sin = sin_exp_func(x_interp, pars);
        y_err = sin_exp_err_model(pars);
        hold on
        x_rad = (x_interp * 180)/pi;
        plot(x_rad, y_sin, lines{mdisp});
        
        [max_y, ind] = max(y_sin);
        [min_y, ind_min] = min(y_sin);
        
        angle_max = x_rad(ind);
        angle_min = x_rad(ind_min);
        
        pref_tilt(1, mdisp) = angle_max;
        %end calculating max for model--------------------------------------------------------
        
        %old max/min determination
        %m_disp_max(mdisp) = spk_max.x;
        %m_disp_min(mdisp) = spk_min.x;
        
        leg{mdisp} = sprintf('%1.3f, %d pts out of range', unique_mean_disp(mdisp), number_disp_outside_range_per_meandisp);
        %save data to calculate the averaged response
        mean_disparity_col = zeros(length(angle_and_response_data), 1) + unique_mean_disp(mdisp);
        all_data(length(all_data)+1:length(all_data)+length(angle_and_response_data),:) = mean_disparity_col;
        angle_data(length(angle_data)+1:length(angle_data)+length(angle_and_response_data),:) = angle_and_response_data;
        
        grad_mean(:,1) = grad_mean(:,1) + px;
        grad_mean(:,2) = grad_mean(:,2) + py;
        
        %angle_string = sprintf(' %1.3f %3.2f %3.2f', unique_mean_disp(mdisp), m_disp_max(mdisp), m_disp_min(mdisp));
        angle_string = sprintf(' %1.3f %3.2f %3.2f', unique_mean_disp(mdisp), angle_max, angle_min);
        angle_out = strcat(angle_out, angle_string);
        
        printcurves = 1;
        if printcurves == 1
            %print out each individual tuning curve for origin
            pathsize = size(PATH,2) - 1;
            while PATH(pathsize) ~='\'	%Analysis directory is one branch below Raw Data Dir
                pathsize = pathsize - 1;
            end   
            PATHOUT = [PATH(1:pathsize) 'Analysis\Grad\'];
            filesize = size(FILE,2) - 1;
            while FILE(filesize) ~='.'
                filesize = filesize - 1;
            end
            FILEOUT = [FILE(1:filesize) 'model_curve'];
            PATHOUT = 'Z:\Users\Jerry\GradAnalysis\';
            fileid = [PATHOUT FILEOUT];
            printflag = 0;
            if (exist(fileid, 'file') == 0)    %file does not yet exist
                printflag = 1;
            end
            proffid = fopen(fileid, 'a');
            if (printflag)
                fprintf(proffid,'X_fit\tY_Fit\tHDisp\tModelResp\tStdErr\tHDisp\tMeasuredResp\tStdErr\n');
                printflag = 0;
            end
            
            for go=1:length(x_rad)
                if go<=length(px)
                    fprintf(proffid,'%6.2f\t%6.2f\t%6.3f\t%6.2f\t%6.2f\t%6.3f\t%6.2f\t%6.2f\n', x_rad(go), y_sin(go), px(go), py(go), perr(go), px_m(go), py_m(go), perr_m(go));
                else
                    fprintf(proffid,'%6.2f\t%6.2f\n', x_rad(go), y_sin(go));
                end
            end
            fclose(proffid);
        end %end printcurves;
        
    end %end mean disparity
    
    %print out saved responses from measured and model for correlation
    %measurement
    print_responses = 1;
    if (print_responses == 1 & printcurves == 1)
        PATHOUT = 'Z:\Users\jerry\GradAnalysis\';
        outfile = [PATHOUT 'model_and_measured_responses_for_' FILE(1:filesize) 'dat'];
        fid = fopen(outfile, 'a');
        for count_mdisp = 1:length(unique_mean_disp)
            for in_mdisp = start_save(count_mdisp):stop_save(count_mdisp)
                fprintf(fid, '%1.2f\t%3.2f\t%6.2f\t%3.2f\t%6.2f\n', unique_mean_disp(count_mdisp), Save_Responses(in_mdisp, 1), Save_Responses(in_mdisp, 2), Save_Responses(in_mdisp, 3), Save_Responses(in_mdisp, 4));
            end
        end
        fclose(fid);
    end
    
    
    
    %calculate TDI from measured data here:
    %readjust mean disparity responses to fall on the same mean
    %then calc avg TDI
    total_mean = mean(TDIdata(:,2));
    for count_meandisp = 1:length(unique_mean_disp)
        disp_mean = mean(TDIdata(start(count_meandisp):stop(count_meandisp),2));
        difference = total_mean - disp_mean;
        TDIdata(start(count_meandisp):stop(count_meandisp),2) = TDIdata(start(count_meandisp):stop(count_meandisp),2) + difference;
    end
    %hold on
    [TDI, var_term] = compute_DDI(TDIdata(:,1)', TDIdata(:,2)');
    
    
    axhandles = get(gca, 'Children');
    numlines = length(unique_mean_disp);
    if legend_on ==1
        leghandles = [];
        legindex = 1;
        for legendloop = numlines:-1:1
            leghandles(legindex) = axhandles((legendloop*3));
            legindex = legindex + 1;
        end
        
        legend(leghandles, leg, 0);
    end
    
    grad_mean(:,1) = grad_mean(:,1)/length(unique_mean_disp);
    grad_mean(:,2) = grad_mean(:,2)/length(unique_mean_disp);
    
    figure(grad_fig)
    hold on;
    [px, py, perr, spk_max, spk_min] = PlotTuningCurve(grad_mean(:,1), grad_mean(:,2), symbols{mdisp+1}, lines{mdisp+1}, 1, 0);
    
    line = sprintf('  %3.2f       %3.2f', spk_max.x, spk_min.x);
    data_string = strcat(data_string, line);
    
    figure(grad_fig);
    height = axis;
    yheight = height(4);
    hold on
    string = sprintf('File: %s', FILE);
    text(height(1), 0.9*yheight, string, 'FontSize', 8);
    
    do_combo = 1;
    if do_combo == 1
        PATHOUT = 'Z:\Users\jerry\GradAnalysis\';
        outfile = [PATHOUT 'pref_tilt_model_vs_pref_tilt_meas_052303.dat'];
        fid = fopen(outfile, 'a');
        if length(unique_mean_disp) > 1
            %go through each combination of mean disparities and print out the relationship between pref tilts
            for j=1:length(unique_mean_disp)
                pref_out = '';
                if p_val(1,j) < 0.05
                    pref_out = sprintf('\t%1.3f\t%3.2f\t%1.3f\t%3.2f\t%1.4f\t%1.4f', unique_mean_disp(j), pref_tilt(1, j), unique_mean_disp(j), pref_tilt_meas(1, j), TDI);
                    line = sprintf('%s', FILE);
                    pref_out = strcat(line, pref_out);
                    fprintf(fid, '%s', [pref_out]);
                    fprintf(fid, '\r');
                end %end if statement
            end %end 1st mdisp search
            fclose(fid);
        end %end check for multiple mdisp
    end %end check if do_combo
    
    
end %end for each magnitude of disparity
%disparities_and_angle

printme = 1;
if (printme==1)
    %pathsize = size(PATH,2) - 1;
    %while PATH(pathsize) ~='\'	%Analysis directory is one branch below Raw Data Dir
    %    pathsize = pathsize - 1;
    %end   
    PATHOUT = 'Z:\Users\jerry\GradAnalysis\';
    
    line = sprintf('%s', FILE);
    data_string = strcat(line, data_string);
    angle_out = strcat(line, angle_out);
    
    %print grad metrics
    outfile = [PATHOUT 'model_predic_052303.dat'];
    fid = fopen(outfile, 'a');
    fprintf(fid, '%s', [data_string]);
    fprintf(fid, '\r\n');
    fclose(fid);
    
    %print mean disp angles
    outfile = [PATHOUT 'Model_Mean_Disp_angles_052303.dat'];
    fid = fopen(outfile, 'a');
    fprintf(fid, '%s', [angle_out]);
    fprintf(fid, '\r\n');
    fclose(fid);        
end
