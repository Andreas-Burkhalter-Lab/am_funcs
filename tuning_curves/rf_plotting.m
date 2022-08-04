%%%% display rf plots
% updated 2018/12/05 on thermaltake

scatter_ctx_center_vs_rf_center = 0; % scatter plot of ctx location vs rf center location
show_individual_rf_maps = 1; % visualize rf and roi location
    show_roi_location = 1; 
    nplots = 49; 
ctx_location_vs_visual_field = 0; 
    
only_plot_significantly_tuned = 1; 
    
tuningdat.tuning.original_ind = [1:height(tuningdat.tuning)]';
if only_plot_significantly_tuned % only plot rois that are significantly spatially tuned
    rows_to_plot = tuningdat.tuning(tuningdat.tuning.rf_sgnf,:);
else
    rows_to_plot = tuningdat.tuning;
end
nrows_to_plot = height(rows_to_plot);

if ctx_location_vs_visual_field
    % visual field rf centers image
    ff = figure('Position',[454 481 1300 420]);
    subplot(1,2,1)
    scatter(rows_to_plot.rf_center_yx(:,2), rows_to_plot.rf_center_yx(:,1),'g','.')
    hold on
    for i = 1:nrows_to_plot
        if rows_to_plot.rf_sgnf(i)
            textcolor = 'r';
        else
            textcolor = 'b';
        end
        origind = rows_to_plot.original_ind(i);
        plottext = text(rows_to_plot.rf_center_yx(i,2), rows_to_plot.rf_center_yx(i,1), num2str(origind), 'Color', textcolor);
    end
    set(gca,'YDir','reverse')
    hold off
    xlabel('nasal --> lateral')
    ylabel('lower --> upper')
    title('visual field')
    %%% ctx plus rois image
    subplot(1,2,2)
    imagesc(tuningdat.meanImage_pre_rotate)
    colormap gray
    hold on
    if strcmp(tuningdat.stimpars.computer,'IMAGING_VR-PC') % de-rotate
        scatter(rows_to_plot.centeryx_pre_rotate(:,2), rows_to_plot.centeryx_pre_rotate(:,1),'g','.')
    elseif strcmp(tuningdat.stimpars.computer,'ANDREWLAB-PC') % no de-rotation
        scatter(rows_to_plot.centeryx(:,2), rows_to_plot.centeryx(:,1),'g','.')
    end
    for i = 1:nrows_to_plot
        if rows_to_plot.rf_sgnf(i)
            textcolor = 'r';
        else
            textcolor = 'b';
        end
        origind = rows_to_plot.original_ind(i);
        if strcmp(tuningdat.stimpars.computer,'IMAGING_VR-PC') % de-rotate
            plottext = text(rows_to_plot.centeryx_pre_rotate(i,2), rows_to_plot.centeryx_pre_rotate(i,1), num2str(origind), 'Color', textcolor);
        elseif strcmp(tuningdat.stimpars.computer,'ANDREWLAB-PC') % no de-rotation
            plottext = text(rows_to_plot.centeryx(i,2), rows_to_plot.centeryx(i,1), num2str(origind), 'Color', textcolor);
        end
    end
    hold off
end

% scatter plot of ctx location vs rf center location
if scatter_ctx_center_vs_rf_center
    figure
    if strcmp(tuningdat.stimpars.computer,'IMAGING_VR-PC') % de-rotate
        scatter(rows_to_plot.centeryx_pre_rotate(:,2),rows_to_plot.rf_center_yx(:,1)) %%% in ctx image, left = ant, up = lateral
        ylabel('visual field: <--lower to upper-->')
        xlabel('ctx location: <--ant to post-->')
        figure
        scatter(rows_to_plot.centeryx_pre_rotate(:,1),rows_to_plot.rf_center_yx(:,2))
        ylabel('visual field: <--nasal to lat-->')
        xlabel('ctx location: <--med to lat-->') %%% assuming presentation into the right visual field
    elseif strcmp(tuningdat.stimpars.computer,'ANDREWLAB-PC') % no de-rotation
        scatter(rows_to_plot.centeryx(:,1),rows_to_plot.rf_center_yx(:,1))  %%% assuming mouse faced straight at screen, no angle, then: in ctx image, left = ant, up = lateral
        ylabel('visual field: <--lower to upper-->')
        xlabel('ctx location: <--ant to post-->')
        [r p] = corrcoef(rows_to_plot.centeryx(:,1),rows_to_plot.rf_center_yx(:,1)); % cor of visual field and ctx coordinates
        cor_r_visupper_brainpost = r(2,1);
        cor_p_visupper_brainpost = p(2,1);
        figure
        scatter(rows_to_plot.centeryx(:,2),rows_to_plot.rf_center_yx(:,2))
        ylabel('visual field: <--nasal to lat-->')
        xlabel('ctx location: <--lat to med-->') %%% assuming presentation into the right visual field
        [r p] = corrcoef(rows_to_plot.centeryx(:,2),rows_to_plot.rf_center_yx(:,2)); % cor of visual field and ctx coordinates
        cor_r_vislat_brainmed = r(2,1);
        cor_p_vislat_brainmed = p(2,1);
        table(cor_r_visupper_brainpost, cor_p_visupper_brainpost, cor_r_vislat_brainmed, cor_p_vislat_brainmed)
    end
end

if show_individual_rf_maps
    figure
    if show_roi_location % show location of roi next to its rf map
        for i = 1:round(nplots/2)
            % resp map
            resp_ind = 2*i-1;
            splotresp = subplot(round(sqrt(nplots)),round(sqrt(nplots)),resp_ind);
            imagesc(rows_to_plot.resp_rf_mean{i})
            colormap(splotresp,'parula')
            colorbar
            title([num2str(rows_to_plot.original_ind(i)), ', anovap ' num2str(rows_to_plot.anovap_rf(i))])
            % roi location
            roi_ind = 2*i;
            splotroi = subplot(round(sqrt(nplots)),round(sqrt(nplots)),roi_ind);
            if strcmp(tuningdat.stimpars.computer,'IMAGING_VR-PC') % de-rotate
                ee = edge(full(rows_to_plot.roi_image_pre_rotate{i}));
            elseif strcmp(tuningdat.stimpars.computer,'ANDREWLAB-PC') % no de-rotation
                ee = edge(full(rows_to_plot.roi_image_prereg{i}));
            end
            [y x] = find(ee);
            roihandle = imagesc(tuningdat.meanImage_pre_rotate);
            colormap(splotroi,'gray')
            hold on
            scatter(x,y,'r','.')
            hold off
        end
    else % just show mean response maps
        for i = 1:nplots
            subplot(round(sqrt(nplots)),round(sqrt(nplots)),i)
            imagesc(rows_to_plot.resp_mean{i})
            colorbar
            title([num2str(rows_to_plot.original_ind(i)), ', anovap ' num2str(rows_to_plot.anovap_rf(i))])
        end
    end
end




