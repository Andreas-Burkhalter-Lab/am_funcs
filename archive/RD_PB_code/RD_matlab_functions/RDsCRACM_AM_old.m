%RD sCRACM AM unused code
% Last modified by AM 2/26/15. 

%% laminar distribution ratios
% % at end of 'm='
%     L1_to_2_5_ratio = L1sum ./ L2_5sum      %% ratio of L1 to L2-5 input for each cell
%     L1_to_2_5_pv = nanmean(L1_to_2_5_ratio(:,1));      %% averages by cell type for plotting (pv cells)
%     L1_to_2_5_pyr = nanmean(L1_to_2_5_ratio(:,2));     %% averages by cell type for plotting (pyr cells)
%     L1_to_2_5_pv_errorbar = std(L1_to_2_5_ratio(:,1) / sqrt(size(L1_to_2_5_ratio,1)));  %% standard error of L1 / L2-5 ratio for PV cells
%     L1_to_2_5_pyr_errorbar = std(L1_to_2_5_ratio(:,2) / sqrt(size(L1_to_2_5_ratio,1)));  %% standard error of L1 / L2-5 ratio for pyramidal cells

% %  plotting
% if plot_layersum
%    figure
%    subplot(2,1,1)       %% subplot placement
%    bar(extractfield(all_cases,'L1_to_2_5_pv')); hold on %% plot L1 : L2-5 ratio for PV cells for all cases
%    errorbar(1:size(all_cases,2),extractfield(all_cases,'L1_to_2_5_pv'), extractfield(all_cases,'L1_to_2_5_pv_errorbar'),'.') %% plot standard error
%    hold off; set(gca,'XTick',1:2*size(all_cases,2),'XTickLabel',CaseList) %% labels
%    title('PV Cells')
%    subplot(2,1,2)       %% subplot placement
%    bar(extractfield(all_cases,'L1_to_2_5_pyr')); hold on %% plot L1 : L2-5 ratio for PV cells for all cases
%    errorbar(1:size(all_cases,2),extractfield(all_cases,'L1_to_2_5_pyr'), extractfield(all_cases,'L1_to_2_5_pyr_errorbar'),'.') %% plot standard error
%    title('Pyramidal Cells')
%    suptitle('L1 : L2-5 Input Ratio')
%    hold off; set(gca,'XTick',1:size(all_cases,2),'XTickLabel',CaseList) %% labels
% end
    
    
%% center vertical distribution on pixel of peak current
%     plot_vertdist_range = [-6:6]; %% maximum distances in pixels from soma to keep in vert_dist (controls y axis range in vertical distribution plot); large values will cause error
% 
% % in 'for c='
%         vertdistmat = zeros(size(plot_vertdist_range,2));  %% erase previous data in 'vertdistmat', include buffer rows
% % in 'for r='
%             % The 'normalize_input' option will affect this section. 
%             % Vertical centering will not be done for L1 and L2-5 sums because raw images are assumed to be aligned with pia. 
%             [vertpeak(r,c) ~] = find(data{r,c}.mean == min(data{r,c}.mean(:))); %% find row index of the peak mean current (soma)
%             vert_dist{r,c} = abs(sum(data{r,c}.mean,2)); % get vertical distribution for this image; remove negative sign
%             L1sum(r,c) = sum(vert_dist{r,c}(L1_pixelrows,:)); % get sum of input into L1 of this cell before centering on soma
%             L2_5sum(r,c) = sum(vert_dist{r,c}(L2_5_pixelrows,:)); % get sum of input into L2-5 of this cell before centering on soma
%             vert_shift(r,c) = round(size(data{r,c}.mean,1) / 2) - vertpeak(r,c);   %% number of pixels that we need to shift 'vert_dist' downwards to center the peak (soma)
%             vert_dist{r,c} = circshift([zeros(round(size(data{r,c}.mean,1) / 2), 1); vert_dist{r,c}; zeros(round(size(data{r,c}.mean,1) / 2), 1)],vert_shift(r,c)); %% buffer, vertically center on peak
%             vert_dist{r,c} = vert_dist{r,c}([-plot_vertdist_range + 2*round(size(data{r,c}.mean,1) / 2)]); %% keep only the specified range of vert_dist
% %  in plotting
%             herrorbar(all_cases(thiscase).group_vertdist(:,thistype),plot_vertdist_range*pixelheight,...    %% horizontally plot vertical distributions with error bars
%                 all_cases(thiscase).vertdist_errorbars(:,thistype),all_cases(thiscase).vertdist_errorbars(:,thistype));

%% centers the image both horizontally and vertically on pixel of peak current
%             [peakind{r,c}(1) peakind{r,c}(2)] = find(data{r,c}.mean == min(data{r,c}.mean(:)));  %% index of the peak mean current
%             hrz_shift(r,c) = round(size(data{r,c}.mean) / 2) - peakind{r,c};   %% number of rows and columns that we need to hrz_shift the image to center to peak
%             mean{r,c} = [zeros(abs(hrz_shift(r,c)(1)),size(data{r,c}.mean,2)); data{r,c}.mean; zeros(abs(hrz_shift(r,c)(1)),size(data{r,c}.mean,2))]; %% add buffer rows on each side to prevent wraparound in cirhrz_shift
%             mean{r,c} = [zeros(size(mean{r,c},1),abs(hrz_shift(r,c)(2))), mean{r,c}, zeros(size(mean{r,c},1),abs(hrz_shift(r,c)(2)))]; %% add buffer columns on each side to prevent wraparound in cirhrz_shift
%             mean{r,c} = circhrz_shift(mean{r,c},hrz_shift(r,c)); %% horizontally hrz_shift enlarged image to center the peak
%             cropcheck = mean{r,c};
%             cropcheck(1+abs(hrz_shift(r,c)(1)):end-abs(hrz_shift(r,c)(1)),1+abs(hrz_shift(r,c)(2)):end-abs(hrz_shift(r,c)(2))) = NaN; %% mask pixels which will not be cropped
%             if any(any(cropcheck))      %% check to see whether nonzero pixels are in the image regions to be cropped
%                 warning('Warning: some nonzero pixels will be cropped while centering on the peak value')
%             end
%             mean{r,c}(1+abs(hrz_shift(r,c)(1)):end-abs(hrz_shift(r,c)(1)),1+abs(hrz_shift(r,c)(2)):end-abs(hrz_shift(r,c)(2))) %% crop peak-centered image to original size

%% Make bar graphs of input to Layer 1 and Layers 2-5 and ratio for all cell types.
% Asssumes that first column is PV cells and second column is pyramidal cells. 
% if plot_layersum
%    figure
%    subplot(2,1,1)       %% subplot placement
%    bar(extractfield(all_cases,'L1_to_2_5_pv')); hold on %% plot L1 : L2-5 ratio for PV cells for all cases
%    errorbar(1:size(all_cases,2),extractfield(all_cases,'L1_to_2_5_pv'), extractfield(all_cases,'L1_to_2_5_pv_errorbar'),'.') %% plot standard error
%    hold off; set(gca,'XTick',1:2*size(all_cases,2),'XTickLabel',CaseList) %% labels
%    title('PV Cells')
%    subplot(2,1,2)       %% subplot placement
%    bar(extractfield(all_cases,'L1_to_2_5_pyr')); hold on %% plot L1 : L2-5 ratio for PV cells for all cases
%    errorbar(1:size(all_cases,2),extractfield(all_cases,'L1_to_2_5_pyr'), extractfield(all_cases,'L1_to_2_5_pyr_errorbar'),'.') %% plot standard error
%    title('Pyramidal Cells')
%    suptitle('L1 : L2-5 Input Ratio')
%    hold off; set(gca,'XTick',1:size(all_cases,2),'XTickLabel',CaseList) %% labels
% end


%% Version of code ~2/23/15 assuming equal number of cells for each cell type.
            %RD sCRACM AM
            % For analyzing RD's sCRACM data, including creating figures similar to those in Petreanu et al. 2009 Figure 2b-d. 
            % Last modified by AM 2/23/15. 

            clear
            close all hidden

            %% Declare cases to include.
            % List (uncommented) in 'CaseList' the pathways/cases to include.

            CaseList={
            %     'PMtoLMorAL_layer23_CRACM'
            %     'R1'
            %    'PMtoV1_layer5_CRACM'
                'V1toPM_L23_sCRACM'
              'V1toPM_L5_sCRACM'
                 'PMtoV1_L23_sCRACM'
            %     'PMtoV1_L5_sCRACM'
            %    'PMtoLM_L23_sCRACM'
            %   'PMtoLM_L5_sCRACM'
            %   'LMtoPM_L23_sCRACM'
            %   'LMtoPM_L5_sCRACM'
                };

            hrz_center = 1;             %% horizontally center all images on peak current pixel before group-averaging
            normto_celltotal = 1;       %% normalize all pixel current values to fraction of total current input into this cell 
            normto_maxpixel = 0;        %% convert all pixel current values to fraction of maximum value within the image


            celltypes = {'PV';'Pyr'};     %% names of cell types; should correspond to columns in CellList
            pixelheight = 75;           %% height of each pixel in um for plotting
            L1_pixelrows = [1:3];          %% rows of pixels from the top of the image representing Layer 1 (vector)
            L2_4_pixelrows = [4:8];       %% rows of pixels from the top of the image representing Layers 2-4 (vector)
            threshold = 0;              %% threshold for classifying cells as receiving or not receiving superthreshold current

            % Warning options: set to 1 to display the warning.
            hemisphere_warn = 0;        %% display warning when hemisphere is not specified
            vertoffset_warn = 0;        %% display warning if no 'vertoffset' value is specified
            crop_warn = 0;              %% display warning when cropping nonzero pixels during horizontal centering

            % Plotting options
            plot_group_images = 1;        % plot group-averaged images (Petreanu et al. 2009 fig. 2b-c)
            plot_vertdist = 1;           % plot group-averaged vertical input distributions (Petreanu et al. 2009 fig. 2d)
            plot_layersum = 0;           % plot groups of total inputs to L1 and to L2-4
            plot_superthresh = 0;        % plot mean number of superthreshold pixels per image within a cell type


            for m=1:size(CaseList,1)        %% for each case/pathway
                CellList = lists_RD1(CaseList{m});      %% load data from each cell from this case
                for c=1:size(CellList,2)    %% for each type of cell (PV and pyramidal)
                    data{1,c} = eval(char(CellList{1,c}));      %% load first cell so we can access 'mean'-field image dimensions
                    cellsmat = zeros([size(data{1,c}.mean),size(CellList,1)]);  %% erase previous data in 'cellsmat;' assumes all images are the same dimensions
                    vertdistmat = zeros(size(data{1,1}.mean,1));  %% erase previous data in 'vertdistmat', include buffer rows
                    for r=1:size(CellList,1)    %% for each recorded cell of this type

                        %% Load matrices of mean current data for each cell and shift vertically according to vertoffset.
                        % load data
                        data{r,c} = eval(char(CellList{r,c}));
                        if strcmp((data{r,c}.hemisphere),'L')==1
                            data{r,c}.mean = fliplr(data{r,c}.mean);
                        elseif hemisphere_warn
                           disp([char(CellList{r,c}) ' does not specify hemisphere, will assume right hemisphere.'])
                        end

                        % Positive values of 'vertoffset' will shift all pixels in the image downwards, negative values shift upwards. 
                        % No warning is provided when cropping nonzero values.   
                        if isfield(data{r,c},'vertoffset')
                            data{r,c}.mean = circshift([zeros(size(data{r,c}.mean)); data{r,c}.mean; zeros(size(data{r,c}.mean))],data{r,c}.vertoffset); %% buffer then shift
                            data{r,c}.mean = data{r,c}.mean(round((1/3)*size(data{r,c}.mean,1))+1:round((2/3)*size(data{r,c}.mean)),:);   %% crop after shifting to restore original image size
                        else
                            data{r,c}.vertoffset = 0;
                            if vertoffset_warn
                                disp([char(CellList{r,c}) ' does not specify vertical alignment offset, assigning offset = 0.'])
                            end
                        end

                        %Changes Inf and -Inf values in maps to NaN
                        data{r,c}.mean(data{r,c}.mean==Inf)=NaN;
                        data{r,c}.mean(data{r,c}.mean==-Inf)=NaN;
                        meanthresholdedmatrix{r,c}=data{r,c}.mean;
                        meanthresholdedmatrix{r,c}(meanthresholdedmatrix{r,c}==0)=NaN;      %% replace all zero values (not meeting threshold) with NaN
                        sumofmeanthresholdedmatrix(r,c) = nansum(meanthresholdedmatrix{r,c}(:));
                        meannotthresholdedmatrix{r,c} = data{r,c}.mean;
                        sumofmeannotthresholdedmatrix(r,c) = sum(data{r,c}.mean(:));

                        %% Normalize and (optionally) horizontally center.
                        % Creating 'groupmean' below requires all images to be the same dimensions.
                        if numel(find([normto_celltotal normto_maxpixel])) > 1
                            error('More than one normalization method should not be selected.')
                        elseif normto_celltotal        %% if the normto_maxpixel option is selected
                            data{r,c}.mean = -(data{r,c}.mean ./ sum(sum(data{r,c}.mean))); %% normalize currents to total current; keep negative sign for use of 'min' below
                        elseif normto_maxpixel        %% if the normto_maxpixel option is selected
                            data{r,c}.mean = -(data{r,c}.mean ./ min(min(data{r,c}.mean))); %% normalize currents to greatest current; keep negative sign for use of 'min' below
                        end

                        % Horizontally center; this section may generate errors if two pixels tie for peak current value.
                        % Horizontally centering only affects the figure displayed by plot_group_images, not other analyses. 
                        if hrz_center
                            [~, peakind{r,c}] = find(data{r,c}.mean == min(data{r,c}.mean(:)));  %% find column index of the peak mean current
                            hrz_shift(r,c) = round(size(data{r,c}.mean,2) / 2) - peakind{r,c};   %% number of columns that we need to shift the image to the right to center the peak
                            peakcentered_mean{r,c} = [zeros(size(data{r,c}.mean,1),abs(hrz_shift(r,c))), data{r,c}.mean, zeros(size(data{r,c}.mean,1),abs(hrz_shift(r,c)))]; %% add buffer columns on each side to prevent wraparound in cirhrz_shift
                            peakcentered_mean{r,c} = abs(circshift(peakcentered_mean{r,c},[0 hrz_shift(r,c)])); %% horizontally shift enlarged image to center the peak; remove negative sign
                            cropcheck = peakcentered_mean{r,c};
                            cropcheck(:,1+abs(hrz_shift(r,c)):end-abs(hrz_shift(r,c))) = NaN; %% mask pixels which will not be cropped
                            if any(any(cropcheck)) && crop_warn      %% check to see whether nonzero pixels are in the image regions to be cropped
                                 warning('for recording "%s," some nonzero pixels will be cropped while horizontally centering on the peak value.',CellList{r,c})
                            end
                            peakcentered_mean{r,c} = peakcentered_mean{r,c}(:,1+abs(hrz_shift(r,c)):end-abs(hrz_shift(r,c))); %% crop peak-centered image to original size
                            cellsmat(:,:,r) = peakcentered_mean{r,c};   %% use the horizontally centered images when creating the group mean below
                        else
                            cellsmat(:,:,r) = data{r,c}.mean;   %% use the uncentered images when creating the group mean below
                        end

                        %% Get vertical distribution, L1 and L2-5 sums, and total superthreshold pixels. 
                        % The 'normalize_input' option will affect this section. 
                        vert_dist{r,c} = abs(sum(data{r,c}.mean,2)); % get vertical distribution for this image; remove negative sign
                        L1sum(r,c) = sum(vert_dist{r,c}(L1_pixelrows,:)); % get sum of input into L1 of this cell before centering on soma
                        L2_4sum(r,c) = sum(vert_dist{r,c}(L2_4_pixelrows,:)); % get sum of input into L2-5 of this cell before centering on soma
                        vertdistmat(:,r) = vert_dist{r,c};     %% put centered vertical distribution into matrix for group averaging below
                        data{r,c}.mean = abs(data{r,c}.mean);    %% remove negative sign
                    end

                %% Get group averages and error bars.
                % Assumes cell type 1 = PV, cell type 2 = pyramidal
                    groupmean{c} = mean(cellsmat,3);    %% create the group mean for this cell type in this case
                    vertdist_errorbars(:,c) = (std(vertdistmat'))' / sqrt(size(vertdistmat,2));    %% find the standard error for each distance from soma for plotting
                    group_vertdist(:,c) = mean(vertdistmat,2);  %% create the group mean centered vertical distribution for this cell type in this case
                    superthresh(r,c) = sum(sum(data{r,c}.mean > threshold));      %% number of pixels with currents exceeding threshold; note +/- sign
                    end
                L1_PV_groupsum = sum(L1sum(:,1)); L2_4_PV_groupsum = sum(L2_4sum(:,1));     %% laminar distribution group sums for PV cells
                L1_pyr_groupsum = sum(L1sum(:,2)); L2_4_pyr_groupsum = sum(L2_4sum(:,2));     %% laminar distribution group sums for PV cells
                L1_PV_errorbars = std(L1sum(:,1)) / sqrt(size(L1sum,1));    %% find standard error of summed L1 input to PV cells
                L1_pyr_errorbars = std(L1sum(:,2)) / sqrt(size(L1sum,1));    %% find standard error of summed L1 input to pyramidal cells
                L2_4_PV_errorbars = std(L2_4sum(:,1)) / sqrt(size(L2_4sum,1));    %% find standard error of summed L2-5 input to PV cells
                L2_4_pyr_errorbars = std(L2_4sum(:,2)) / sqrt(size(L2_4sum,1));    %% find standard error of summed L2-5 input to pyramidal cells
                superthresh_PV = mean(superthresh(:,1)); superthresh_pyr = mean(superthresh(:,2));   %% mean superthreshold by cell type for plotting
                superthresh_errorbars_PV = std(superthresh(:,1)) / sqrt(size(superthresh,1));  %% standard error of number of superthreshold pixels for PV cells
                superthresh_errorbars_pyr = std(superthresh(:,2)) / sqrt(size(superthresh,1));  %% standard error of number of superthreshold pixels for pyramidal cells

                casename = CaseList(m); %% create case name variable for v2struct
                clear vars; vars = who; vars = ['fieldNames'; vars(cellfun(@isempty,strfind(vars,'all_cases')))]; %% names of variables to put into all_cases
                all_cases(m) = v2struct(vars); 
            end 

            %% Create figures
            set(0,'DefaultTextInterpreter','none'); %% suptitle doesn't take an 'Interpreter' argument so it must be changed here.

            % Create group images - fig. 2b/c from Petreanu et al. 2009 
            if plot_group_images
                figure 
                for thiscase = 1:size(CaseList,1)
                    for thistype = 1:size(celltypes,1)    %% for each cell type
                        subplot(size(celltypes,1),size(CaseList,1),thiscase+size(CaseList,1)*(thistype-1))       %% subplot placement
                        imagesc(all_cases(thiscase).groupmean{thistype}); title(sprintf('%s - %s',CaseList{thiscase},celltypes{thistype}),'FontSize',8);
            %             xlabel(sprintf('Pixel Length = %g um',pixelheight)); ylabel(sprintf('Pixel Length = %g um',pixelheight))
                    end
                end
                    if hrz_center
                            suptitle('Group Means (Horizontally Centered)');
                        else
                            suptitle('Group Means (Not Centered'); %% title
                    end
                set(gcf,'Position',[0 300 620 500])    %% figure position and size 
                set(gca,'Fontsize',13);  %% last number controls font size
            end

            % Create vertical distribution plot - fig. 2d (red line) from Petreanu et al. 2009
            if plot_vertdist
                figure
                for thiscase = 1:size(CaseList,1)   % one column for each case
                    for thistype = 1:size(celltypes,1)    %% one subplot for each cell type
                        subplot(size(celltypes,1),size(CaseList,1),thiscase+size(CaseList,1)*(thistype-1))       %% subplot placement
                        herrorbar(all_cases(thiscase).group_vertdist(:,thistype),fliplr(1:size(all_cases(thiscase).group_vertdist(:,thistype))),...%% horizontally plot vertical distributions with error bars
                            all_cases(thiscase).vertdist_errorbars(:,thistype),all_cases(thiscase).vertdist_errorbars(:,thistype));
                        xlimits = xlim; xlim([0 xlimits(2)]);   %% set x-axis minimum to zero
                        title(sprintf('%s - %s',CaseList{thiscase},celltypes{thistype}),'FontSize',8); %% last value controls font size
                        xlabel('Normalized Input'); %%ylabel('Distance From Soma (um)')
                    end
                end
                suptitle('Input Vertical Distribution'); set(gcf,'Position',[650 300 620 500])    %% title, position, window size
            end

            % Make bar graphs of input to Layer 1 and Layers 2-5 for all cell types.
            % Asssumes that first column is PV cells and second column is pyramidal cells. 
            % 'Barweb' may not work on Matlab 2014b; use "hold on; errorbar(x,y,error,'rx')" instead. 
            if plot_layersum
               % PV Cells
               figure
               subplot(2,1,1)       %% subplot placement
               barweb([(extractfield(all_cases,'L1_PV_groupsum'))',(extractfield(all_cases,'L2_4_PV_groupsum'))'],...
                   [(extractfield(all_cases,'L1_PV_errorbars'))',(extractfield(all_cases,'L2_4_PV_errorbars'))'],...
                   [],CaseList','PV Cells Input Received in L1 and L2-4',[],'Normalized Input',[],[],{'L1','L2-4'}); %% plot PV L1 and L2-4 sums
               set(gca,'Fontsize',9)

               % Pyramidal Cells
               subplot(2,1,2)       %% subplot placement
               barweb([(extractfield(all_cases,'L1_pyr_groupsum'))',(extractfield(all_cases,'L2_4_pyr_groupsum'))'],...
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

            set(0,'DefaultTextInterpreter','Tex');  %% reset DefaultTextInterpreter
            all_cases = rmfield(all_cases,{'hemisphere_warn','crop_warn','cropcheck','cellsmat','vertdistmat','m','r','c',...
                'plot_group_images','plot_vertdist','plot_layersum'}); %% clear unneeded variables

            % clearvars -except all_cases;        %% clear workspace except for structure containing separate analyses of cases



            % fix superthresh
            % n = number of cell
            % option to eliminate repeats

            %% To do:
            % [idea is that PV cells get similar input from larger area]
            % check barweb error bars - too small? 
            % image not imagesc
            % show what gets cropped
            % ?interpolate to create hi-res image?