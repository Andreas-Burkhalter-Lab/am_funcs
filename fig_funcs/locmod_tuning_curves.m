 %%%% plot tuning curves of cells in different quantiles to show different level of locomotion modulation
% filetable used in tuning_stats must be loaded
 % updated 2020/8/18 on thermaltake
 
 %%%%% NOTE FOR IMPORTING INTO ILLUSTRATOR: open the .svg or .eps exported from Matlab as its own file, then copy into the larger figure
 %%%%% ....... inserting the vector image with the 'place' command may result in distorting of the tuning curve lines
 
%  close all
 
 exemplar_cells_excel = 'F:\thesis\exemplar_cells.xlsx';
%  excel_rows = [2 10 17]; % rows to plot
 excel_rows = [2 10 17 24 28 36 47 54 70]; % rows to plot
override_locomotion_thresh = 0; % use loc thresh specified in this script, rather than in excel
    loc_thresh_vals = [0, logspace(-5,-3,8)]; % plot all of these log thresh values
    subplot_row_col_per_cell = [3 3];
show_plots = 1;  
     axes_line_width = 4.5;
     curve_line_width = 4; 
     use_subplots = 1; %%% plot all cells in subplots in a single figure
        subplot_row_col = [3 3]; 
        subplot_spacing_yx = [0.1 0.08]; 
    axes_numbers_bold = 'bold'; %%% 'bold' or 'normal'
    axis_font_size = 15;  
    xticklabels_linear = 1; % force xticks to be of form 0.1 rather than 10^-1 (does not effect tick location or axis scale)
    sf_xticks = 10.^[-2 -1 0];
    tf_xticks = [0.2 1 10];
    orient_xticks = [0 100 200 300]; 
    fig_outerposition = [0.05 0.05 0.6 0.9]; 
    %     fig_outerposition = [0.05 0.05 0.3 0.4]; % [left bottom width height]
    show_subtitle = 0; 
    show_suptitle = 0; 
    font_name = 'Arial';

    save_image = 1; % save the final figure as an .svg
        savename = 'F:\thesis\figs\fig 2 -- tuning curves.svg'; 
store_tuningdat = 1; % load all tuningdat structs into workspace

    
if ~exist('filetable','var')
    filetable = excel2table('F:\analysis_master.xlsx'); % load master file list
    filetable = filetable(filetable.analyze_plane == 1,:); % keep only planes that have be marked for analysis
    filetable.analyze_plane = []; 
end
keyvars = {'centeryx_prereg','quantile'}; % vars from excel that will also be found in tuningdat.tuning
frontvars = {'pakan_loc_index','sf_anovap','tf_anovap','orient_anovap','locm_cor_r'}; % variables to bring to front of table
 celltable = readtable(exemplar_cells_excel);
 celltable = celltable(excel_rows-1, :); % keep only specified cells
 celltable.pakan_loc_index = []; %%% make join.m usable
 ncells = height(celltable);
 celltable.centeryx_prereg(:,2) = celltable.centeryx_prereg2; % format
 celltable.centeryx_prereg2 = [];
 if override_locomotion_thresh; nplots_per_cell=subplot_row_col_per_cell(1)*subplot_row_col_per_cell(2); else nplots_per_cell=1;  end
 
 tuning_plot_pars.LineWidth = curve_line_width;
 tuning_plot_pars.show_title = show_subtitle;

 
 % plot tuning curve for cell selected cells
 clear temptable
 for icell = 1:ncells 
     % load tuningdat if it's not already loaded
     store_tuningdat_completed = vardefault('store_tuningdat_completed',false); 
     tuningdat_stored = vardefault('tuningdat_stored',table); 
     if ~store_tuningdat_completed
         if ~exist('filetable_row','var') || ~[filetable_row == find( filetable.sub == celltable.subject(icell) & strcmp(filetable.day,char(celltable.day(icell))) & filetable.plane == celltable.plane(icell))]
             filetable_row = find( filetable.sub == celltable.subject(icell) & strcmp(filetable.day,char(celltable.day(icell))) & filetable.plane == celltable.plane(icell)); 
             load([filetable.directory{filetable_row}, filesep, filetable.tuning_file{filetable_row}])
         end
         matchrow = find( tuningdat.tuning.centeryx_prereg(:,1) == celltable.centeryx_prereg(icell,1) &... % find the specified cell from this session
                    tuningdat.tuning.centeryx_prereg(:,2) == celltable.centeryx_prereg(icell,2) );
         if store_tuningdat
            tuningdat_stored.tuningdat{icell} = tuningdat;
            tuningdat_stored.matchrow(icell) = matchrow; 
         end
         if ncells < 2
            temptable(icell,:) = join(celltable(icell,:),tuningdat.tuning,'Keys',keyvars); % ?causes errors when run with multiple cells
         end
     end
     if show_plots
         if icell == 1 || ~[use_subplots || override_locomotion_thresh]
             tuningfig = figure('units','normalized','outerposition',fig_outerposition); 
         end
         if use_subplots % put all cells in one figure
            subtightplot(subplot_row_col(1), subplot_row_col(2), icell, subplot_spacing_yx)
         end
         for iplot = 1:nplots_per_cell
             if nplots_per_cell > 1
                 subplot(subplot_row_col_per_cell(1), subplot_row_col_per_cell(2), iplot)
             end
             if override_locomotion_thresh
                 this_loc_thresh = loc_thresh_vals(iplot); 
             else
                 this_loc_thresh = celltable.locomotion_thresh(icell);
             end
             tuning_plot_pars = vardefault('tuning_plot_pars',struct);
             if store_tuningdat_completed
                plot_tuning_curve(tuningdat_stored.tuningdat{icell}, tuningdat_stored.matchrow(icell), celltable.param{icell}, this_loc_thresh, tuning_plot_pars); 
             elseif ~store_tuningdat_completed
                plot_tuning_curve(tuningdat, matchrow, celltable.param{icell}, this_loc_thresh, tuning_plot_pars); 
             end
             set(gca,'FontName', font_name);
            if show_subtitle
                title(['loc_thresh ', num2str(this_loc_thresh)])
            end
            xlabel('')
            ylabel('')
            hLeg=findobj(gcf,'type','legend');
                set(hLeg,'visible','off') % turn off legend
            hax = get(gca); 

            set(gca,'Box','off')
            set(gca,'linewidth',axes_line_width)
            set(gca,'FontSize',axis_font_size)
            set(gca,'FontWeight',axes_numbers_bold)
            switch celltable.param{icell}
                case 'sf'
                    hax.XAxis.TickValues = sf_xticks;
                case 'tf'
                    hax.XAxis.TickValues = tf_xticks;
                case 'orient'
                    hax.XAxis.TickValues = orient_xticks;
            end
            set(gca,'LooseInset',get(gca,'TightInset')+[0 0 0.005 0]) % crop borders
            if xticklabels_linear
                hax.XAxis.TickLabels = cellstr(num2str(hax.XAxis.TickValues')); % force x axis labels out of exponent mode
            end
         end
         if show_suptitle
            suptitle([num2str(celltable.subject(icell)), ' ', char(celltable.day(icell))])
         end
     end
 end
try celltable = temptable; clear temptable; end
try celltable = movevars(celltable,frontvars,'After','locomotion_thresh'); end
if store_tuningdat
    store_tuningdat_completed = 1; 
end
if save_image
    if exist(savename,'file')
        delete(savename)
    end
    saveas(tuningfig, savename, 'svg')
end