%%% plot multiple response parameters of an roi from gcamp recording
% function plot_roi_response(tuningdata, tablerow)
%
%%% updated 19/3/4 on thermaltake

function plot_roi_response(tuningdata, tablerow)

separate_high_locomotion_trials = 1; 

%% tuning curve plotting
if ~isempty(tuningdata.tuningpars.stimpars.sf_vals)
    subplot(2,3,1)
    plot_tuning_curve(tuningdata, tablerow, 'sf', separate_high_locomotion_trials)
end

if ~isempty(tuningdata.tuningpars.stimpars.tf_vals)
    subplot(2,3,2)
    plot_tuning_curve(tuningdata, tablerow, 'tf', separate_high_locomotion_trials)
end
    
if ~isempty(tuningdata.tuningpars.stimpars.orient_vals)
    subplot(2,3,3)
    plot_tuning_curve(tuningdata, tablerow, 'orient', separate_high_locomotion_trials)
end
    
%% rf map
splotresp = subplot(2,3,4);
rf_on_stimscreen = NaN(tuningdata.rfdat.stimpars.commonvars.screenstats.height, tuningdata.rfdat.stimpars.commonvars.screenstats.width); 
rf_stim_radius_pix = tuningdata.rfdat.stimpars.diam_pix / 2; 
% get the coords of corner pixels covered by the rf grid
grid_min_y = round(tuningdata.rfdat.stimpars.stim_centers_yx_pix{1,1}(1) - rf_stim_radius_pix);
grid_max_y = round(tuningdata.rfdat.stimpars.stim_centers_yx_pix{end,end}(1) + rf_stim_radius_pix);
grid_min_x = round(tuningdata.rfdat.stimpars.stim_centers_yx_pix{1,1}(2) - rf_stim_radius_pix);
grid_max_x = round(tuningdata.rfdat.stimpars.stim_centers_yx_pix{end,end}(2) + rf_stim_radius_pix);
y_pix_per_row = round([grid_max_y - grid_min_y] / tuningdata.rfdat.stimpars.rows);
x_pix_per_column = round([grid_max_x - grid_min_x] / tuningdata.rfdat.stimpars.columns);
%%% grid responses projected on the screen image (one element per screen pixel)
resp_rf_mean_upsampled = repelem(tuningdata.tuning.resp_rf_mean{tablerow}, y_pix_per_row, x_pix_per_column); 
rf_on_stimscreen(grid_min_y:grid_min_y+size(resp_rf_mean_upsampled,1)-1, grid_min_x:grid_min_x+size(resp_rf_mean_upsampled,2)-1) = ...
    resp_rf_mean_upsampled; %%% add response map to screen image

% reconstruct warped coordinates
eye2screen_center = sqrt( tuningdata.stimpars.eye2screen_top_bottom^2 - (tuningdata.stimpars.screen_height/2)^2 ); %because eye2screen_edge^2=eye2screen_center^2 + (screen_length/2)^2
eye2screen_center_pix = eye2screen_center * (tuningdata.stimpars.screenstats.height / tuningdata.stimpars.screen_height); % convert to pix from inches
scrnperpx = round(tuningdata.stimpars.screenstats.width/2); % perpendicular to screen through eye must pass through this x-value; usually screen center
scrnperpy = round(tuningdata.stimpars.screenstats.height/2); % perpendicular to screen through eye must pass through this y-value; usually screen center
perp2stimcenterx = scrnperpx - tuningdata.stimpars.stim_center_yx(2); % x pixels between screen center and stim center
perp2stimcentery = scrnperpy - tuningdata.stimpars.stim_center_yx(1); % y pixels between screen center and stim center
[xcentermesh, ycentermesh] = meshgrid(1:tuningdata.stimpars.screenstats.width,1:tuningdata.stimpars.screenstats.height);
perp2meshx = scrnperpx-xcentermesh;
perp2meshy = scrnperpy-ycentermesh;
apt_fullscreen = deg2rad(tuningdata.stimpars.diam_minmax(2)/2) > acos( (eye2screen_center_pix^2 + perp2stimcenterx*perp2meshx + perp2stimcentery*perp2meshy) ./...
sqrt( (eye2screen_center_pix^2+perp2stimcenterx^2+perp2stimcentery^2)*(eye2screen_center_pix^2+perp2meshx.^2+perp2meshy.^2)) );
aperture_edge_img = edge(apt_fullscreen);
[apt_edge_y, apt_edge_x] = find(aperture_edge_img);

% plot rf against stim location
imagesc(rf_on_stimscreen) %% plot the rf map
colormap(splotresp,'parula')
colorbar
hold on
scatter(apt_edge_x,apt_edge_y,'.','r') %%% show stimulus outline
title(['anovap ' num2str(tuningdata.tuning.anovap_rf(tablerow))])
hold off

%% roi location
splotroi = subplot(2,3,5);
if strcmp(tuningdata.stimpars.computer,'IMAGING_VR-PC') % de-rotate
    ee = edge(full(tuningdata.tuning.roi_image_pre_rotate{tablerow}));
elseif strcmp(tuningdata.stimpars.computer,'ANDREWLAB-PC') % no de-rotation
    ee = edge(full(tuningdata.tuning.roi_image_prereg{tablerow}));
end
[y x] = find(ee);
roihandle = imagesc(tuningdata.meanImage_pre_rotate);
colormap(splotroi,'gray')
hold on
scatter(x,y,'r','.')
hold off

%% dff timecourse
splot_dff = subplot(2,3,6);
load(tuningdata.dff_file);
plot(dff(:,tablerow))
title('dff')