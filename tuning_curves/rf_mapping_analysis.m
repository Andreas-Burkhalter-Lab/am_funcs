%%%%% perform receptive field mapping on all rois by fitting a 2d gaussian
% updated 2020/1/16 on thermaltake

show_rf_plots = 0; 

% for each roi, make table of F response for each trial sorted by stim parameter and trial 
tuning = table(NaN(nrois,2),  false(nrois,1),NaN(nrois,1),NaN(nrois,2),NaN(nrois,2),       NaN(nrois,2),  cell(nrois,1), cell(nrois,1),      cell(nrois,1),       NaN(nrois,2),           cell(nrois,1),...
    'VariableNames',{'centeryx','rf_sgnf', 'anovap_rf','rf_center_yx','rf_center_yx_pix', 'rf_size_yx_deg', 'resp_rf',   'resp_rf_mean',      'roi_image_prereg','centeryx_pre_rotate','roi_image_pre_rotate'});
respmatrix = NaN(stimpars.rows,stimpars.columns,stimpars.repetitions); % record all responses for a given roi in 3d matrix
% for 2d gauss fitting
[meanspikes_grid_x meanspikes_grid_y] = meshgrid(1:stimpars.columns, 1:stimpars.rows); 
[x_predict, y_predict] = meshgrid(... %% make an x-y plane for projecting the receptive field onto the stimulus screen, by interpolating and expanding grid index values
      [ ([1:stimpars.commonvars.screenstats.width] - stimpars.stim_centers_yx_pix{1}(2) )./stimpars.h_spacing_pix ] + 1,...  %%% x values; add 1 so that stim grid coord [1,1] has value 1
      [ ([1:stimpars.commonvars.screenstats.height] -stimpars.stim_centers_yx_pix{1}(1) )./stimpars.v_spacing_pix ] + 1 ); %%% y values

tuning.rf_fit = cell(nrois,1);
tuning.rf_surf = cell(nrois,1);
leftmost_x_pix = min(stimpar_sets.stim_center_x);
topmost_y_pix = min(stimpar_sets.stim_center_y);
stimpars.deg_per_pix = stimpars.h_spacing / stimpars.h_spacing_pix; % value used in stim presentation, only true for closest point on screen to eye

wbar = waitbar(0,['Computing receptive field maps...']);
for iroi = 1:nrois 
    record_roi_location()
    if isvalid(wbar)
        waitbar(iroi/nrois,wbar)
    end
    % tabulate responses
    tuning.resp_rf{iroi} = respmatrix;
    for irow = 1:stimpars.rows
        for icolumn = 1:stimpars.columns
            match = find(stimpar_sets.row == irow & stimpar_sets.column == icolumn);
            resp_this_rowcol = stimpar_sets.dff_during_stim(match,iroi); % all responses to this row col combination
            tuning.resp_rf{iroi}(irow,icolumn,1:length(resp_this_rowcol)) = reshape(resp_this_rowcol,1,1,[]); % rotate into third dimension, put into resp matrix
        end
    end
    tuning.resp_rf_mean{iroi} = nanmean(tuning.resp_rf{iroi},3); % take means of stim at the same locations
    % significance test - simple anova using grid locations as groups, not using spatial information
    table_for_anova = [reshape(tuning.resp_rf{iroi}, stimpars.rows*stimpars.columns,[])]';
    tuning.anovap_rf(iroi) = anova1(table_for_anova,[],'off');
    tuning.rf_sgnf(iroi) = tuning.anovap_rf(iroi) < 0.05;
    
    % best-responding grid location
    maxval = max(tuning.resp_rf_mean{iroi}(:));
    [tuning.prefindyx(iroi,1), tuning.prefindyx(iroi,2)] = find(tuning.resp_rf_mean{iroi} == maxval);
    
    if pars.fit_tuning_functions
        % fit 2d gaussian
        tuning.rf_fit{iroi} = struct;
        [tuning.rf_fit{iroi}.fitresult, tuning.rf_fit{iroi}.zfit, tuning.rf_fit{iroi}.fiterr, tuning.rf_fit{iroi}.zerr,...
            tuning.rf_fit{iroi}.resnorm, tuning.rf_fit{iroi}.rr] = ...  % fit the grid response data with a 2D Gaussian
            fmgaussfit(meanspikes_grid_x, meanspikes_grid_y, tuning.resp_rf_mean{iroi});
        tuning.rf_center_yx(iroi,:) = fliplr(tuning.rf_fit{iroi}.fitresult(5:6));
        % Generate a 3D surface from the 2D Gaussian fitted parameters for plotting.
        rf_surf = single( tuning.rf_fit{iroi}.fitresult(7) + tuning.rf_fit{iroi}.fitresult(1)*exp(... % single to save space, may want to compress more or not save at all
            -(((x_predict-tuning.rf_fit{iroi}.fitresult(5)).*cosd(tuning.rf_fit{iroi}.fitresult(2))+...
            (y_predict-tuning.rf_fit{iroi}.fitresult(6)).*sind(tuning.rf_fit{iroi}.fitresult(2)))./tuning.rf_fit{iroi}.fitresult(3)).^2-...
            ((-(x_predict-tuning.rf_fit{iroi}.fitresult(5)).*sind(tuning.rf_fit{iroi}.fitresult(2))+...
            (y_predict-tuning.rf_fit{iroi}.fitresult(6)).*cosd(tuning.rf_fit{iroi}.fitresult(2)))./tuning.rf_fit{iroi}.fitresult(4)).^2) );
        % find pixels on the stim screen covering the cell's receptive field, defined as those within 2 standard deviations of 2d Gaussian peak (from Gao et al. 2010)
        standdev_x = tuning.rf_fit{iroi}.fitresult(3) / sqrt(2); % in units of grid spaces - standdev of ROTATED x-axis of the curve, not x-axis of the stim screen
        %%% plug in [rf peak plus 2 standard deviations] (x and y) to the 2d gaussian function to find the threshold response for points within the receptive field
        % compute threshold on the rotated curve by setting theta to zero; y coordinate set to peak y
        tuning.rf_fit{iroi}.rf_threshold_response = tuning.rf_fit{iroi}.fitresult(7) + tuning.rf_fit{iroi}.fitresult(1)*exp(...
            -[((2*standdev_x))./tuning.rf_fit{iroi}.fitresult(3)].^2);
        tuning.rf_fit{iroi}.rf_image = rf_surf > tuning.rf_fit{iroi}.rf_threshold_response; 
%         tuning.rf_surf{iroi} = rf_surf; %%% generally don't save rf surf, because it multiplies file size by about 10x
        % get the width and height of the receptive field - using screen x and y directions, will NOT match the sizes of the rotated ellipse used to draw the RF image
        [rf_ycoord, rf_xcoord] = find(tuning.rf_fit{iroi}.rf_image); % rf_coord will be empty if rf is projected to fall completely off screen
        if ~isempty(rf_ycoord) && min(rf_ycoord) > 1 && max(rf_ycoord) < size(tuning.rf_fit{iroi}.rf_image, 1) % if rf doesn't go past the top or bottom of the screen, else leave as nan
            tuning.rf_size_yx_deg(iroi,1) = stimpars.deg_per_pix * range(rf_ycoord); % get rf height
        end
        if ~isempty(rf_xcoord) && min(rf_xcoord) > 1 && max(rf_xcoord) < size(tuning.rf_fit{iroi}.rf_image, 2) % if rf doesn't go past the left or right edges of the screen, else leave as nan
            tuning.rf_size_yx_deg(iroi,2) = stimpars.deg_per_pix * range(rf_ycoord); % get rf width
        end
        clear rf_surf
    end
end
if isvalid(wbar)
    close(wbar)
end   

tuning.rf_size_yx_deg(tuning.rf_size_yx_deg==0) = NaN; %% make sure that rf size was not incorrectly assigned as zero
tuning.rf_center_yx_pix(:,1) = stimpars.v_spacing_pix * [tuning.rf_center_yx(:,1) - 1] + topmost_y_pix; % convert grid location y center to pix
tuning.rf_center_yx_pix(:,2) = stimpars.h_spacing_pix * [tuning.rf_center_yx(:,2) - 1] + leftmost_x_pix; % convert grid location x center to pix
res.tuning = tuning;

if show_rf_plots
    rf_plotting();
end

