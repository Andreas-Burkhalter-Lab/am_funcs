%RD sCRACM AM
% For analyzing RD's sCRACM data, including creating figures similar to those in Petreanu et al. 2009 Figure 2b-d.
% Assumes that first cell type is PV, second is pyramidal.
% Last modified by AM 8/19/15 on recording comp. 

clear
close all hidden
global group_images_scalebar_labels cbar % for fixing the plot outside of this function

%% Declare cases to include.
% List (uncommented) in 'CaseList' the pathways/cases to include.

CaseList={
%   'PMtoLMorAL_layer23_CRACM'
%   'R1'
%   'PMtoV1_layer5_CRACM'
%   'V1toPM_L23_sCRACM'
%   'V1toPM_L5_sCRACM'
%   'PMtoV1_L23_sCRACM'
%   'PMtoV1_L5_sCRACM'
%   'PMtoLM_L23_sCRACM'
   'PMtoLM_L5_sCRACM'
%   'LMtoPM_L23_sCRACM'
%   'LMtoPM_L5_sCRACM'
    };

hrz_center = 1;             %% horizontally center all images on peak current pixel before group-averaging
normto_pairtotal = 0;       %% normalize pixel current values to fraction of total current input avereged across the pair
normto_celltotal = 0;       %% normalize all pixel current values to fraction of total current input into this neuron 
normto_maxpixel = 0;        %% convert all pixel current values to fraction of maximum value within the image
normto_pairmax = 0;         %% normalize pixel values to maximum pixel value  
no_repeats = 1;            %% do not include repeat copies of neurons (only one cell per neuron name)
vertdist_subjects_are = 'neurons';   %% choose 'neurons' or 'pixels'; unit to treat as a subject for creating vertdist error bars

celltypes = {'PV';'Pyr'};     %% names of cell types; should correspond to columns in CellList
pixelheight = 75;           %% height of each pixel in um for plotting
L1_pixelrows = [1:3];          %% rows of pixels from the top of the image representing Layer 1 (vector)
L2_4_pixelrows = [4:16];       %% rows of pixels from the top of the image representing Layers 2-4 (vector)
L4_pixelrows = [3:5];
threshold = 0;              %% threshold for classifying cells as receiving or not receiving superthreshold current
manova_on = 0;         %% set to 1 to perform multivariate anova on input vertical distrbition between cases
spanova_on = 0;         %% set to 1 to perform split-plot anova on input vertical distrbition between cases 
anova_rows = [3:8];     %% anova_rows determines pixel rows to use as variables for both manova and split-plot anova (spanova)

% Warning options: set to 1 to display the warning.
hemisphere_warn = 0;        %% display warning when hemisphere is not specified
vertoffset_warn = 0;        %% display warning if no 'vertoffset' value is specified
crop_warn = 0;              %% display warning when cropping nonzero pixels during horizontal centering

% Plotting options
plot_group_images={1,'scaled',12}  ; %1st arg: plot group mean; 2nd arg: set to 'scaled' to scale colors within each plot; 3rd arg: resolution enhancement factor (=1 for original image
     group_images_scalebar_labels = [0 0.25 0.5 0.75 1]; % linearly-spaced labels for scale bar; can be any number of labels
plot_vertdist = {1, 'overlay', 'nogrid'};  % argument 1 = 1 for plotting on; argument 2 = 'overlay' to plot all cases on the same graph; argument 3 = 'grid' to display grid lines
plot_vertdist_difference = 0;       % plot difference between input into each row in first pathway vs. second pathway
plot_layersum = 0;           % plot groups of total inputs to L1 and to L2-4
plot_superthresh = 1;        % plot mean number of superthreshold pixels per image within a cell type
plot_3rows = {0, [1 2 3]};           %% first argument enables plotting of 3 rows of image as 3-D coordinate; secdon argument controls which rows are plotted (must be 3 integers)
plot_L1_rows_sums = 0;         % plot mean input into each row of layer 1, treating each row of layer 1 as a subject

for m=1:size(CaseList,1)        %% for each case/pathway
    CellList = lists_RD1(CaseList{m});      %% load data from each cell from this case
    %% Turn columns of CellList into cell arrays to allow for different subject numbers. 
    for c = 1:size(CellList,2)   %% for each type of cell (PV and pyramidal)
        tempcelllist{c}{1} = CellList{1,c};     %% first row cannot be a repeat of higher rows
        data{c}{1} = eval(char(CellList{1,c})); %% first row cannot be a repeat of higher rows
        for r = 2:size(CellList,1)    %% for each cell listed (including repeats)
            if ~(no_repeats && any(strcmp(CellList{r,c},CellList(1:r-1,c))));   %% as long as no_repeats is off or the name isn't a repeat of a preceiding row...
                data{c}{size(tempcelllist{c},2)+1} = eval(char(CellList{r,c})); %% ...include this cell data as next on the list
                tempcelllist{c}{size(tempcelllist{c},2)+1} = CellList{r,c};     %% ...make this cell name the next name on the list
            end
        end
    end
    CellList = tempcelllist; clear tempcelllist;    %% remove repeats from CellList if applicable
    origdata = data; % store an original copy so that we can change the values of data
    %% Analysis
    for c=1:size(CellList,2)    %% for each type of cell (PV and pyramidal)
        cellsmat = zeros([size(data{c}{1}.mean),size(CellList{c},2)]);  %% erase previous data in 'cellsmat;' assumes all images are the same dimensions
        L1_nrows = length(L1_pixelrows);
        L1_rows_sums{c} = NaN(L1_pixelrows* size(CellList{c},2),1); % initialize
        for r=1:size(CellList{c},2)    %% for each recorded cell of this type
            % Hemisphere check
            if strcmp((data{c}{r}.hemisphere),'L')==1
                data{c}{r}.mean = fliplr(data{c}{r}.mean);
            elseif hemisphere_warn
               disp([char(CellList{c}{r}) ' does not specify hemisphere, will assume right hemisphere.'])
            end
            
            % Positive values of 'vertoffset' will shift all pixels in the image downwards, negative values shift upwards. 
            % No warning is provided when cropping nonzero values.   
            if isfield(data{c}{r},'vertoffset')
                data{c}{r}.mean = circshift([zeros(size(data{c}{r}.mean)); data{c}{r}.mean; zeros(size(data{c}{r}.mean))],data{c}{r}.vertoffset); %% buffer then shift
                data{c}{r}.mean = data{c}{r}.mean(round((1/3)*size(data{c}{r}.mean,1))+1:round((2/3)*size(data{c}{r}.mean)),:);   %% crop after shifting to restore original image size
            else
                data{c}{r}.vertoffset = 0;
                if vertoffset_warn
                    disp([char(CellList{c}{r}) ' does not specify vertical alignment offset, assigning offset = 0.'])
                end
            end
            %Changes Inf and -Inf values in maps to NaN
            data{c}{r}.mean(data{c}{r}.mean==Inf)=NaN;
            data{c}{r}.mean(data{c}{r}.mean==-Inf)=NaN;
        
            % Normalize and (optionally) horizontally center. ('unique' is an alternative.)
            % Creating 'groupmean' below requires all images to be the same dimensions.
            if numel(find([normto_pairtotal normto_celltotal normto_maxpixel normto_pairmax])) > 1
                error('More than one normalization method should not be selected.')
            elseif normto_pairtotal 
                if any(size(CellList{1}) ~= size(CellList{2}))
                    error('Number of PV and Pyr cells must be equal to normalize using pairs. (''no_repeats'' might be active.)')
                else
                    pair_avg_celltotal = mean([sum(sum(origdata{1}{r}.mean)), sum(sum(origdata{2}{r}.mean))]); % mean total current into the pair
                    data{c}{r}.mean = -(data{c}{r}.mean ./ pair_avg_celltotal); % normalize to pair average total input
                end
            elseif normto_celltotal        %% if the normto_maxpixel option is selected
                data{c}{r}.mean = -(data{c}{r}.mean ./ sum(sum(data{c}{r}.mean))); %% normalize currents to total current; keep negative sign for use of 'min' below
            elseif normto_maxpixel        %% if the normto_maxpixel option is selected
                data{c}{r}.mean = -(data{c}{r}.mean ./ min(min(data{c}{r}.mean))); %% normalize currents to greatest current; keep negative sign for use of 'min' below
            elseif normto_pairmax
                if any(size(CellList{1}) ~= size(CellList{2}))
                    error('Number of PV and Pyr cells must be equal to normalize using pairs. (''no_repeats'' might be active.)')
                else
                pairmax = max([min(min(origdata{1}{r}.mean)), min(min(origdata{2}{r}.mean))]);
                data{c}{r}.mean = -(data{c}{r}.mean ./ pairmax);
                end
            end
            
            % Horizontally center; this section may generate errors if two pixels tie for peak current value.
            % Horizontally centering only affects the figure displayed by plot_group_images, not other analyses. 
            if hrz_center
                [~, peakind{c}{r}] = find(data{c}{r}.mean == min(data{c}{r}.mean(:)));  %% find column index of the peak mean current
                hrz_shift{c}{r} = round(size(data{c}{r}.mean,2) / 2) - peakind{c}{r};   %% number of columns that we need to shift the image to the right to center the peak
                peakcentered_mean{c}{r} = [zeros(size(data{c}{r}.mean,1),abs(hrz_shift{c}{r})), data{c}{r}.mean, zeros(size(data{c}{r}.mean,1),abs(hrz_shift{c}{r}))]; %% buffer columns
                peakcentered_mean{c}{r} = abs(circshift(peakcentered_mean{c}{r},[0 hrz_shift{c}{r}])); %% horizontally shift enlarged image to center the peak; remove negative sign
                cropcheck = peakcentered_mean{c}{r};
                cropcheck(:,1+abs(hrz_shift{c}{r}):end-abs(hrz_shift{c}{r})) = NaN; %% mask pixels which will not be cropped
                if any(any(cropcheck)) && crop_warn      %% check to see whether nonzero pixels are in the image regions to be cropped
                     warning('for recording "%s," some nonzero pixels will be cropped while horizontally centering on the peak value.',CellList{c}{r})
                end
                peakcentered_mean{c}{r} = peakcentered_mean{c}{r}(:,1+abs(hrz_shift{c}{r}):end-abs(hrz_shift{c}{r})); %% crop peak-centered image to original size
                cellsmat(:,:,r) = peakcentered_mean{c}{r};   %% use the horizontally centered images when creating the group mean below
            else
                cellsmat(:,:,r) = data{c}{r}.mean;   %% use the uncentered images when creating the group mean below
            end
            
            % Get vertical distribution, L1 and L2-5 sums, L1 rows sums, and total superthreshold pixels. 
            % The 'normalize_input' option will affect this section. 
            vert_dist{c}(:,r) = abs(sum(data{c}{r}.mean,2)); % get vertical distribution for this image; remove negative sign
            superthresh_L1{c}(r) = numel(find(data{c}{r}.mean(L1_pixelrows,:)));      %% number of superthreshold pixels in L1 of this neuron
            superthresh_L2_4{c}(r) = numel(find(data{c}{r}.mean(L2_4_pixelrows,:)));      %% number of superthreshold pixels in L2-4 of this neuron
            data{c}{r}.mean = abs(data{c}{r}.mean);    %% remove negative sign
            L1_rows_sums{c}(L1_nrows*(r-1)+1 : L1_nrows*r,1) = vert_dist{c}(L1_pixelrows,r); % one sum per L1 row of this neuron
            
        end
        % Get group averages and error bars.
        % Assumes cell type 1 = PV, cell type 2 = pyramidal
        groupmean{c} = mean(cellsmat,3);    %% create the group-mean image for this cell type in this case
        hrz_concat{c} = reshape(cellsmat,size(cellsmat,1),[],1);    %% concatenate all cell images horizontally to remove 3rd dimension
        L1_pix_all{c} = hrz_concat{c}(L1_pixelrows,:);     %% concatenate all L1 pixels from cells of this type
        L1_pix_nonzero{c} = L1_pix_all{c}(L1_pix_all{c} > 0);       %% concatenate all nonzero L1 pixels
        L2_4_pix_all{c} = hrz_concat{c}(L2_4_pixelrows,:);     %% concatenate all L2-4 pixels from cells of this type 
        L2_4_pix_nonzero{c} = L2_4_pix_all{c}(L2_4_pix_all{c} > 0);       %% concatenate all nonzero L2-4 pixels
        group_vertdist(:,c) = mean(vert_dist{c},2);  %% create the group-mean vertical distribution for this cell type in this case
        L1sum{c} = sum(vert_dist{c}(L1_pixelrows,:)); % get sums of input into L1 of cells of this type
        L2_4sum{c} = sum(vert_dist{c}(L2_4_pixelrows,:)); % get sums of input into L2-4 of cells of this type
        L4_pix_all{c} = hrz_concat{c}(L4_pixelrows,:);
        L4_pix_nonzero{c} = L4_pix_all{c}(L4_pix_all{c} > 0);  
        
%         if ~ischar(vertdist_subjects_are) || ~any(strcmp(vertdist_subjects_are,{'neurons','rows','pixels',}));
%             vertdist_subjects_are = questdlg('What units will be treated as subjects for creating error bars?',...
%                 'Subject Units','neurons','rows','pixels','neurons');
%         end
        switch vertdist_subjects_are
            case 'neurons'
                vertdist_errorbars(:,c) = (std(vert_dist{c}'))' / sqrt(size(vert_dist{c},2)); % SEM for each row, treating each neuron as a subject
%             case 'row'
%                 
            case 'pixels'
                vertdist_errorbars(:,c) = (std(hrz_concat{c}'))' / sqrt(size(hrz_concat{c},2)); % SEM for each row, treating each pixel as a subject
        end       
        superthresh{c} = superthresh_L1{c} + superthresh_L2_4{c};      %% number of pixels with currents exceeding threshold; note +/- sign
    end
    L1_PV_groupmean = mean(L1sum{1}); L2_4_PV_groupmean = mean(L2_4sum{1});     %% laminar distribution group means for PV cells
    L1_pyr_groupmean = mean(L1sum{2}); L2_4_pyr_groupmean = mean(L2_4sum{2});     %% laminar distribution group means for pyramidal cells
    L1_PV_errorbars = std(L1sum{1}) / sqrt(size(L1sum{1},2));    %% find standard error of summed L1 input to PV cells
    L1_pyr_errorbars = std(L1sum{2}) / sqrt(size(L1sum{2},2));    %% find standard error of summed L1 input to pyramidal cells
    L2_4_PV_errorbars = std(L2_4sum{1}) / sqrt(size(L2_4sum{1},2));    %% find standard error of summed L2-4 input to PV cells
    L2_4_pyr_errorbars = std(L2_4sum{2}) / sqrt(size(L2_4sum{2},2));    %% find standard error of summed L2-4 input to pyramidal cells
    superthresh_PV = mean(superthresh{1}); superthresh_pyr = mean(superthresh{2});   %% mean superthreshold by cell type for plotting
    superthresh_errorbars_PV = std(superthresh{1}) / sqrt(size(superthresh{1},2));  %% standard error of number of superthreshold pixels for PV cells
    superthresh_errorbars_pyr = std(superthresh{2}) / sqrt(size(superthresh{2},2));  %% standard error of number of superthreshold pixels for pyramidal cells
    L1_rows_sums_PV = mean(L1_rows_sums{1}); % mean pixel sum per L1 row for PV neurons
    L1_rows_sums_pyr = mean(L1_rows_sums{2}); % mean pixel sum per L1 row for pyramidal neurons
    L1_rows_sums_errorbars_PV = std(L1_rows_sums{1}) / sqrt(size(L1_rows_sums{1},1)); % standard error of L1 rows sums for PV neurons
    L1_rows_sums_errorbars_pyr = std(L1_rows_sums{2}) / sqrt(size(L1_rows_sums{2},1)); % standard error of L1 rows sums for PV neurons
    
    casename = CaseList(m); %% create case name variable for v2struct
    clear vars; vars = who; vars = ['fieldNames'; vars(cellfun(@isempty,strfind(vars,'all_cases')))]; %% names of variables to put into all_cases
    all_cases(m) = v2struct(vars);  %% put analysis results into the next element of all_cases
    clear data vert_dist superthresh_L1 superthresh_L2_4 groupmean hrz_concat L1_pix_all L2_4_pix_all L1_pix_nonzero L2_4_pix_nonzero group_vertdist L1sum L2_4sum...
        vertdist errorbars superthresh L1_PV_errorbars L1_pyr_errorbars L2_4_PV_errorbars L2_4_pyr_errorbars superthresh_errorbars_PV superthresh_errorbars_pyr
end  

% Create 'padded' pixel vectors by the concatenating nonzero-pixel vectors 
% with zeros so that they are all as large as the largest one among all 
% cases within a given cell type and layergroup (L1 vs L2-4).
L1maxpix(1) = max(cellfun(@(x) size(x{1},1),extractfield(all_cases,'L1_pix_nonzero')));    %% max size of L1_pix_nonzero among all cases for PV cells
L1maxpix(2) = max(cellfun(@(x) size(x{2},1),extractfield(all_cases,'L1_pix_nonzero')));    %% max size of L1_pix_nonzero among all cases for pyramidal cells
L2_4maxpix(1) = max(cellfun(@(x) size(x{1},1),extractfield(all_cases,'L2_4_pix_nonzero')));    %% max size of L2_4_pix_nonzero among all cases for PV cells
L2_4maxpix(2) = max(cellfun(@(x) size(x{2},1),extractfield(all_cases,'L2_4_pix_nonzero')));    %% max size of L2_4_pix_nonzero among all cases for pyramidal cells
for m = 1:size(all_cases,2)   %% for each case
    for c = 1:size(all_cases(1,1).celltypes,1)  %% for each cell type
        all_cases(m).L1_pix_padded{c} = [all_cases(m).L1_pix_nonzero{c}; zeros(L1maxpix(c)-size(all_cases(m).L1_pix_nonzero{c},1),1)];  %% create padded L1 pixel vector
        all_cases(m).L2_4_pix_padded{c} = [all_cases(m).L2_4_pix_nonzero{c}; zeros(L2_4maxpix(c)-size(all_cases(m).L2_4_pix_nonzero{c},1),1)];  %% create padded L2_4 pixel vector        
    end
end

%% Create figures
set(0,'DefaultTextInterpreter','none'); %% suptitle doesn't take an 'Interpreter' argument so it must be changed here.

if plot_group_images{1}
    figure 
    cbar = NaN(length(celltypes),length(CaseList));
    maxpix = extractfield(all_cases,'groupmean'); maxpix = max(max(cell2mat(cellfun(@cell2mat,maxpix,'UniformOutput',false)))); %% maximum pixel value for heat map scaling
    for thiscase = 1:size(CaseList,1)
        for thistype = 1:size(celltypes,1)    %% for each cell type
            subplot(size(celltypes,1),size(CaseList,1),thiscase+size(CaseList,1)*(thistype-1))       %% subplot placement
            if strcmp(plot_group_images{2},'scaled')
                imagesc(imresize(all_cases(thiscase).groupmean{thistype},plot_group_images{3})); 
                title(sprintf('%s - %s',CaseList{thiscase},celltypes{thistype}),'FontSize',8);                
            else
                image(imresize(64/maxpix*all_cases(thiscase).groupmean{thistype},plot_group_images{3})); 
                title(sprintf('%s - %s',CaseList{thiscase},celltypes{thistype}),'FontSize',8); 
            end
%             xlabel(sprintf('Pixel Length = %g um',pixelheight)); ylabel(sprintf('Pixel Length = %g um',pixelheight))            
            set(gca,'Fontsize',10);  %% last number controls font size
            cbar(thistype,thiscase) = colorbar; set(gca,'pos',get(gca,'pos')+[0.01 0 0 0]); %% units will be normalized to 64
            set(cbar(thistype,thiscase),'YTick',linspace(0,64,length(group_images_scalebar_labels)));
            set(cbar(thistype,thiscase),'YTickLabel',group_images_scalebar_labels);
        end
    end
    if hrz_center; suptitle('Group Means (Horizontally Centered)');
    else suptitle('Group Means (Not Centered'); %% title
    end
    set(gcf,'Position',[0 300 620 500])    %% figure position and size 
end


% Create vertical distribution plot - fig. 2d (red line) from Petreanu et al. 2009
if plot_vertdist{1}
    figure;     
    set(gcf,'Position',[650 300 620 500]); %suptitle('Input Vertical Distribution');     %% title, position, window size
    for thiscase = 1:size(CaseList,1)   % one column for each case
        for thistype = 1:size(celltypes,1)    %% one subplot for each cell type
            if strcmp(plot_vertdist{2},'overlay')
                vdist_forplotting(:,thiscase,thistype) = all_cases(1,thiscase).group_vertdist(:,thistype);
                verror_forplotting(:,thiscase,thistype) = all_cases(1,thiscase).vertdist_errorbars(:,thistype);
                subplot(size(celltypes,1),1,thistype)       %% overlay plots from different cases
                if thiscase == size(CaseList,1);    %% if this is the final case
                    herrorbar(vdist_forplotting(:,:,thistype),repmat((fliplr(1:size(vdist_forplotting,1)))',1,size(vdist_forplotting,2)),...
                        verror_forplotting(:,:,thistype),verror_forplotting(:,:,thistype))
                    title(sprintf(celltypes{thistype}),'FontSize',9); %% last value controls font size
                    legend(CaseList,'Interpreter','none','FontSize',8);  %% last value controls font size
                end    %% 
            else    %% if each case is to be plotted on a separate graph
                subplot(size(celltypes,1),size(CaseList,1),thiscase+size(CaseList,1)*(thistype-1))       %% new column for each case\
                herrorbar(all_cases(thiscase).group_vertdist(:,thistype),fliplr(1:size(all_cases(thiscase).group_vertdist(:,thistype))),...%% horizontally plot vertical distributions with error bars
                    fliplr(all_cases(thiscase).vertdist_errorbars(:,thistype)),fliplr(all_cases(thiscase).vertdist_errorbars(:,thistype)));
                title(sprintf('%s - %s',CaseList{thiscase},celltypes{thistype}),'FontSize',8); %% last value controls font size
            end
            if strcmp(plot_vertdist{3},'grid'); grid; end
            xlimits = xlim; xlim([0 xlimits(2)]);   %% set x-axis minimum to zero
            xlabel('Normalized Input'); %%ylabel('Distance From wm (um)')
        end
    end
end

% Plot differences between vertical distributions in first pathway vs. second pathway. 
% Only designed to work with 2 pathways. Error bars will be added if each
% pathway has an equal number of cells of a given type. 
if plot_vertdist_difference
    figure;      
    for thistype = 1:size(celltypes,1)    %% one subplot for each cell type
        for thiscase = 1:size(CaseList,1)   % one column for each case
            vdist_forplotting(:,thiscase,thistype) = all_cases(1,thiscase).group_vertdist(:,thistype);           
            if length(all_cases(1,1).CellList{thistype}) == length(all_cases(1,2).CellList{thistype})   %% if each case has an equal number of neurons of this type
                vdist_allsubs(:,:,thiscase,thistype) = all_cases(1,thiscase).vert_dist(1,thistype); 
            end
            if thiscase == size(CaseList,1);    %% if this is the final case
                vdist_dif_group(:,thistype) = vdist_forplotting(:,2,thistype) - vdist_forplotting(:,1,thistype);  %% difference between vertdist between the 2 pathways
                subplot(size(celltypes,1),1,thistype)      
                hold all
                vdist_dif_group_pos = vdist_dif_group(:,thistype); vdist_dif_group_pos( vdist_dif_group_pos < 0 ) = 0; % keep positive values
                vdist_dif_group_neg = vdist_dif_group(:,thistype); vdist_dif_group_neg( vdist_dif_group_pos > 0 ) = 0; % keep negative values
                barh(flipud(vdist_dif_group_pos), 'r')  % rows where case 2 > case 1
                barh(flipud(vdist_dif_group_neg), 'b')  % rows where case 2 > case 1
                title(sprintf(celltypes{thistype}),'FontSize',9); %% last value controls font size
                xlabel(sprintf('%s Input - %s Input (Normalized Units)', CaseList{2}, CaseList{1})); %%ylabel('Distance From Soma (um)')
                   set(gca,'Fontsize',9); set(gcf,'Position',[1300 300 620 500])    %% font size, position, window size 
                if length(all_cases(1,1).CellList{thistype}) == length(all_cases(1,2).CellList{thistype})   %% if each case has an equal number of neurons of this type
                    vdist_dif_allsubs{thistype} = vdist_allsubs(:,:,2,thistype) - vdist_allsubs(:,:,1,thistype);
                    vdist_dif_errorbars(thistype) = std(vdist_dif_allsubs{thistype} / sqrt(size(vdist_allsubs,2))); %% size(vdist_allsubs,2) equals number of neurons of this type
                    herrorbar(vdist_dif_group(:,thistype), 1:length(vdist_dif_group(:,thistype)), vdist_dif_errorbars(thistype), vdist_dif_errorbars(thistype))
                end
            end
        end
    end
end
                        
% Make bar graphs of input to Layer 1 and Layers 2-5 for all cell types.
% Asssumes that first column is PV cells and second column is pyramidal cells. 
% 'Barweb' may not work on Matlab 2014b; use "hold on; errorbar(x,y,error,'rx')" instead. 
if plot_layersum
   % PV Cells
   figure
   subplot(2,1,1)       %% subplot placement
   barweb([(extractfield(all_cases,'L1_PV_groupmean'))',(extractfield(all_cases,'L2_4_PV_groupmean'))'],...
       [(extractfield(all_cases,'L1_PV_errorbars'))',(extractfield(all_cases,'L2_4_PV_errorbars'))'],...
       [],CaseList','PV Cells Input Received in L1 and L2-4',[],'Normalized Input',[],[],{'L1','L2-4'}); %% plot PV L1 and L2-4 sums
   set(gca,'Fontsize',9)
   
   % Pyramidal Cells
   subplot(2,1,2)       %% subplot placement
   barweb([(extractfield(all_cases,'L1_pyr_groupmean'))',(extractfield(all_cases,'L2_4_pyr_groupmean'))'],...
       [(extractfield(all_cases,'L1_pyr_errorbars'))',(extractfield(all_cases,'L2_4_pyr_errorbars'))'],...
       [],CaseList','Pyramidal Cells Input Received in L1 and L2-4',[],'Normalized Input',[],[],{'L1','L2-4'}); %% plot PV L1 and L2-4 sums
   set(gca,'Fontsize',9); set(gcf,'Position',[1300 300 620 500])    %% font size, position, window size 
end

% Make bar graphs of total number of superthreshold pixels all cell types.
% Asssumes that first column is PV cells and second column is pyramidal cells. 
if plot_superthresh
   % PV Cells
   figure
   subplot(2,1,1)       %% subplot placement - PV cells
   bar(extractfield(all_cases,'superthresh_PV')); hold on %% plot mean superthreshold pixels for PV cells for all cases
   errorbar(1:size(all_cases,2),extractfield(all_cases,'superthresh_PV'), extractfield(all_cases,'superthresh_errorbars_PV'),'.') %% plot standard error
   title('PV Cells','Fontsize',12)  %% last number controls font size of title
   set(gca,'XTick',1:size(all_cases,2),'XTickLabel',CaseList) %% case name labels
   set(gca, 'Fontsize',10) %% font size of case names
   
   % Pyramidal Cells
   subplot(2,1,2)       %% subplot placement - pyramidal cells
   bar(extractfield(all_cases,'superthresh_pyr')); hold on %% plot mean superthreshold pixels for PV cells for all cases
   errorbar(1:size(all_cases,2),extractfield(all_cases,'superthresh_pyr'), extractfield(all_cases,'superthresh_errorbars_pyr'),'.') %% plot standard error
   title('Pyramidal Cells','Fontsize',12)  %% last number controls font size of titles
   set(gca,'XTick',1:size(all_cases,2),'XTickLabel',CaseList) %% case name labels
   set(gca, 'Fontsize',10) %% font size of case names
   hold off;suptitle('Number of Pixels Receiving Input')
end

% Plot mean values of inputs from 3 rows as 3-D coordinates, comparing
% across both cases and neuron types. 
if plot_3rows{1}
    figure, hold all
    for thiscase = 1:size(CaseList,1)   % for each case
        for thistype = 1:size(celltypes,1)    %% for each neuron type type
            plot3(all_cases(1,thiscase).group_vertdist(plot_3rows{2}(1),thistype),all_cases(1,thiscase).group_vertdist(plot_3rows{2}(2),thistype),...
                all_cases(1,thiscase).group_vertdist(plot_3rows{2}(3),thistype),'o','DisplayName',sprintf('%s - %s',CaseList{thiscase},celltypes{thistype}));
        end
    end
    lgnd = legend('show'); set(lgnd,'Interpreter','none'); grid on;
    xlabel(sprintf('Mean Input Into Row %g',plot_3rows{2}(1))),ylabel(sprintf('Mean Input Into Row %g',plot_3rows{2}(2))),...
        zlabel(sprintf('Mean Input Into Row %g',plot_3rows{2}(3)))
end

% Make bar graphs of pixel value sums per row of layer 1, treating each row
if plot_L1_rows_sums
   % PV Cells
   figure
   subplot(2,1,1)       %% subplot placement - PV cells
   bar(extractfield(all_cases,'L1_rows_sums_PV')); hold on %% plot mean superthreshold pixels for PV cells for all cases
   errorbar(1:size(all_cases,2),extractfield(all_cases,'L1_rows_sums_PV'), extractfield(all_cases,'L1_rows_sums_errorbars_PV'),'.') %% plot standard error
   title('PV Neurons','Fontsize',12)  %% last number controls font size of title
   set(gca,'XTick',1:size(all_cases,2),'XTickLabel',CaseList) %% case name labels
   set(gca, 'Fontsize',10) %% font size of case names
   
   % Pyramidal Cells
   subplot(2,1,2)       %% subplot placement - pyramidal cells
   bar(extractfield(all_cases,'L1_rows_sums_pyr')); hold on %% plot mean superthreshold pixels for PV cells for all cases
   errorbar(1:size(all_cases,2),extractfield(all_cases,'L1_rows_sums_pyr'), extractfield(all_cases,'L1_rows_sums_errorbars_pyr'),'.') %% plot standard error
   title('Pyramidal Neurons','Fontsize',12)  %% last number controls font size of titles
   set(gca,'XTick',1:size(all_cases,2),'XTickLabel',CaseList) %% case name labels
   set(gca, 'Fontsize',10) %% font size of case names
   
   set(gcf,'Position',[1300 300 620 500]); % right side
   hold off;suptitle('Mean Layer 1 Row Sums (Rows as Subjects)')
end

set(0,'DefaultTextInterpreter','Tex');  %% reset DefaultTextInterpreter
all_cases = rmfield(all_cases,{'hemisphere_warn','crop_warn','vertoffset_warn','cropcheck','cellsmat','m','r','c',...
    'plot_group_images','plot_vertdist','plot_layersum',}); %% clear unneeded variables
    

%% MANOVA and Split-Plot Anova
% Compare vertical input distributions between cases.
% In the MANOVA, each row of pixels is treated as a separate variable. 
if manova_on
    manova_vdist = cell(1,size(celltypes,1)); manova_grouplabels = cell(1,size(celltypes,1));       %% create blank so we can add to these variables later
    for thistype = 1:size(celltypes,1)    %% separate neuron types into different cells
        for thiscase = 1:size(CaseList,1)   % vertically concatenate vert_dist from all cases into manova_vdist
            manova_vdist{thistype} = [manova_vdist{thistype}; (all_cases(1,thiscase).vert_dist{thistype}(anova_rows,:))'];
            manova_grouplabels{thistype} = [manova_grouplabels{thistype}; repmat(all_cases(1,thiscase).casename,size(all_cases(1,thiscase).vert_dist{thistype},2),1)];
        end
        [manova_dim(thistype) manova_p(thistype) manova_stats{thistype}] = manova1(manova_vdist{thistype},manova_grouplabels{thistype});
    end
    manova_p        %% display manova p-values
end

if spanova_on 
    spanova_vdist = cell(1,size(celltypes,1)); spanova_groupingvars = cell(1,size(celltypes,1));  %% create blank so we can add to these variables later
    NestingStructure = [ 0 0 0; 0 0 0; 1 0 0];      %% Grouping Variable Nesting Hierarchy - from Mike Rieger's code
    ModelStructure = [1 0 0; 0 0 1; 1 1 0; 0 1 0];  %% from Mike Rieger's code
    for thistype = 1:size(celltypes,1)
        spanova_groupingvars{thistype} = cell(1,3);    %% for each neuron type, one cell for each of the grouping variables CASE, PIXELROW, and NEURON ID
        for thiscase = 1:size(CaseList,1)        
            spanova_vdist{thistype} = [spanova_vdist{thistype}; reshape(all_cases(1,thiscase).vert_dist{thistype}(anova_rows,:),[],1)]; %vcat all vertdists from all neurons of this type
            spanova_groupingvars{thistype}{1,1} = [spanova_groupingvars{thistype}{1};...
                repmat(all_cases(1,thiscase).casename,length(anova_rows)* size(all_cases(1,thiscase).vert_dist{thistype},2),1)]; %% case names
            spanova_groupingvars{thistype}{1,2} = [nominal(spanova_groupingvars{thistype}{2});...
                nominal(repmat((1:length(anova_rows))',size(all_cases(1,thiscase).vert_dist{thistype},2),1))]; %% pixel rows
            spanova_groupingvars{thistype}{1,3} = [ nominal(spanova_groupingvars{thistype}{3});...
                nominal((reshape(repmat(all_cases(1,thiscase).CellList{thistype},length(anova_rows),1),1,[]))') ];     %% neuron IDs
        end
        [spanova_p(:,thistype) spanova_table{thistype}] = anovan(spanova_vdist{thistype},spanova_groupingvars{thistype},'random',3,... %calculate split-plot ANOVA - from Mike Rieger's code
                        'nested',NestingStructure,'model',...
                        ModelStructure,'varnames',{'Group','Distance','NeuronID'},'display','off');
        spanova_table{thistype} = spanova_table{thistype}(:,[1,2,3,5,11,10,6,7]);  %% clean up split-plot ANOVA results - from Mike Rieger's code       
    end
    spanova_p = spanova_p(1,:)
end

clearvars -except all_cases celltypes CaseList manova_dim manova_p manova_stats spanova_p spanova_table;   %% clear workspace except for analyses of cases and manova/spanova results

% plot showing differences between pathways by layer

%% To do:
% look at multcompare

