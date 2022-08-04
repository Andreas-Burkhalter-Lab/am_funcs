%%% laminar_density_plotting: plot optical density laminar distribution from LA, dlGN, LP to POR from multiple cases
%
% interpolate values for missing sections, get means and standard error for each layer and pathway
%
% updated 2020/2/21 on thermaltake 


filelist = 'C:\Users\Burkhalter Lab\Documents\anatomy_paper\lam_dens_file_list.xlsx';
axis_font_size = 20; 
axes_line_width = 2;
axes_numbers_bold = 'bold'; %%% 'bold' or 'normal'
legend_labels = {'LA-->POR'; 'dLGN-->POR'; 'LP-->POR'};
line_width = 4;
box_on_off = 'off';
section_thickness_um = 40; % for converting section number to depth
horizontal_plot = 1;
fig_width = 1000;% pixels
fig_height = 950; % pixels
ylimits = [0 600];
xticks_hrzplot = linspace(0, 1, 5); % x ticks for horizontal plot
xticklabels_hrzplot = {'0','','0.5','','1'};
subplot_extra_spacing = 0.06; %%% space to add between subplots
save_lam_density_fig = 0; 
    savename_lam_density = 'fig_lam_density';

colororder = [... % plot line colors
0    0    1;...
0.4660    0.6740    0.1880;...
1   0    0;...
0.4940    0.1840    0.5560;...
0.3010    0.7450    0.9330;...
0.6350    0.0780    0.1840;...
0.9290    0.6940    0.1250;...
];

sec_to_depth = @(sec,thick)thick*sec - thick;

results_by_case = areaDensity(filelist); %%% analyze and tabulate laminar density data

%% sort results by pathway and subject
unqpaths = unique(results_by_case.pathway,'stable');
npaths = length(unqpaths);
lamdata = table(unqpaths,cell(npaths,1),cell(npaths,1), 'VariableNames', ...
                {'pathway', 'subjects', 'data'});
for ipath = 1:npaths
    thispath = unqpaths{ipath};
    lamdata.pathway{ipath} = thispath; 
    this_path_rows = strcmp(results_by_case.pathway,thispath);
    thesesubs = unique(results_by_case.sub(this_path_rows));
    n_subs_this_path = length(thesesubs);
    lamdata.subjects{ipath} = thesesubs; 
    n_secs_this_path = max(results_by_case.section(this_path_rows)); % highest number section among cases in this pathway
    nancol = NaN(n_secs_this_path,1);
    nanmat = NaN(n_secs_this_path,n_subs_this_path); 
    lamdata.data{ipath} = table(nancol, nanmat,    nanmat,        nanmat,         nancol,     nancol,   'VariableNames', ...
                            {'section','depth_um','intens','intens_pre_interp','mean_intens','sem_intens'});
    lamdata.data{ipath}.section = [1:n_secs_this_path]';
    for isub = 1:n_subs_this_path
        this_sub_rows = this_path_rows & results_by_case.sub == thesesubs(isub); % rows match this subject and path
        this_sub_sections = results_by_case.section(this_sub_rows); 
        this_sub_normedIntens = results_by_case.normedIntens(this_sub_rows); 
        lamdata.data{ipath}.intens_pre_interp(this_sub_sections,isub) = this_sub_normedIntens; %%% add data to the table for this subjects' intensities
        lamdata.data{ipath}.intens(:,isub) = interp1(this_sub_sections, this_sub_normedIntens, lamdata.data{ipath}.section, 'linear'); % interpolate missing section intensities
        lamdata.data{ipath}.mean_intens = mean(lamdata.data{ipath}.intens,2);
        lamdata.data{ipath}.sem_intens = std(lamdata.data{ipath}.intens,0,2) / sqrt(n_subs_this_path);
    end
    lamdata.data{ipath}.depth_um = sec_to_depth(lamdata.data{ipath}.section, section_thickness_um); % assign approximate depths
end

 %% make plot           
fig_lam_density = figure;
hold on
for ipath = 1:npaths
    hsubplot(ipath) = subplot(1,npaths,ipath);
    if horizontal_plot
        hplot = herrorbar( lamdata.data{ipath}.mean_intens, lamdata.data{ipath}.depth_um, lamdata.data{ipath}.sem_intens);%
        hplot(1).Color = colororder(ipath,:);
        hplot(2).Color = colororder(ipath,:);
%         hplot(1).Color = [0 .6 0];
%         hplot(2).Color = [0 .6 0];
        hplot(1).LineWidth = line_width;
        hplot(2).LineWidth = line_width;
        plot_formatting()
        set(gca,'YDir','reverse')
%         xlabel('Optical density')
        set(gca,'XTick',xticks_hrzplot)
        set(gca,'XTickLabels',xticklabels_hrzplot)
%         ylabel('Depth (µm)')
    elseif ~horizontal_plot
        errorbar(lamdata.data{ipath}.depth_um, lamdata.data{ipath}.mean_intens, lamdata.data{ipath}.sem_intens, 'Color',colororder(ipath,:),'LineWidth',line_width)
        %     ylabel('Optical density')
        xlabel('Depth (µm)')
    end

    set(gca,'defaultAxesColorOrder',colororder)


end

% add spacing betewen subplots
hsubplot(1).Position(1) = hsubplot(1).Position(1) - subplot_extra_spacing;
hsubplot(3).Position(1) = hsubplot(3).Position(1) + subplot_extra_spacing;

% legend(legend_labels)

set(fig_lam_density,'Renderer', 'painters', 'Position', [50 50 fig_width fig_height]) % set figure length and width
if save_lam_density_fig
    print(fig_lam_density,savename_lam_density, '-dtiffn','-r300') %%% save image as file
end
