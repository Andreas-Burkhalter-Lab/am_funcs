%%%% find threshold intensity for patches
%%% updated 18/8/20 on thermaltake

%%% blurring
H = fspecial('disk',diskblurradius_pix); %%% create filter; arguments set blur radius
imroi = im;
imroi(~roi) = NaN;
patchdata.imroi  = imroi; 
if pars.blur_using_only_roi
    imroi = double(imroi);  % must be double to have nans
    imroiblurred = nanconv(imroi,H,'nanout','edge');
    imblurred = [];
else
    imblurred = imfilter(im,H,'replicate'); %% 
    imroiblurred = imblurred;
    imroiblurred(~roi) = NaN;
end
imroivals = imroiblurred(roi);

minInRoi = min(imroivals);
maxInRoi = max(imroivals);

switch pars.threshMode
    case 'intensityQuantiles' % get threshold for each quantile
        nQuantiles = pars.nQuantiles;
        roisorted = sort(imroivals);
        quantvals = quantile(imroivals,1/nQuantiles:1/nQuantiles:1);
        quantile_table = table(NaN(nQuantiles,1),cell(nQuantiles,1),cell(nQuantiles,1),cell(nQuantiles,1),'Variablenames',{'pixIntensThresh','quantimage','edges','pixyx'});
        quant_levels_img = zeros(size(im));
        thresh =  quantvals(nQuantiles-pars.topNQuantilesArePatches);
        % divide image into intensity quantiles
        for iquant = 1:nQuantiles
            if iquant == 1
                quantile_table.quantimage{iquant} = imroiblurred <= quantvals(iquant);
            else % eliminate overlaps between quantiles
                quantile_table.quantimage{iquant} = imroiblurred <= quantvals(iquant) & imroiblurred > quantvals(iquant-1);
                quantile_table.pixIntensThresh(iquant) = quantvals(iquant-1); %% minimum pix intens to equal or exceed to be included in this quant group
            end
            [quantile_table.pixyx{iquant}(:,1), quantile_table.pixyx{iquant}(:,2)] = find(quantile_table.quantimage{iquant});
            quant_levels_img = quant_levels_img + iquant*double(quantile_table.quantimage{iquant});
            quantile_table.edges{iquant} = edge(quantile_table.quantimage{iquant}); % get edges of this quantile
        end
        patchdata.quantile_table = quantile_table;
        patchdata.quant_levels_img = quant_levels_img;  
    case 'fractionOfPopulation'
        roisorted = sort(imroivals);
        threshpoint = max([1, length(roisorted) - round(fracPixSuperthresh * length(roisorted))]);
        thresh = roisorted(threshpoint); 
        patchdata.threshpoint = threshpoint;
        patchdata.threshFraction_fractionOfMax = double(thresh)/double(max(roisorted)); % convert for comparison with other threshmode
        threshinc = maxInRoi * pars.raisethresh_increment;
    case 'fractionOfMaxMinusMin'
        thresh =  minInRoi + [maxInRoi-minInRoi] * threshFraction;
        threshinc = [maxInRoi-minInRoi] * pars.raisethresh_increment;
    case 'fractionOfMax'
        thresh =  maxInRoi * threshFraction;
        threshinc = maxInRoi * pars.raisethresh_increment;
end
    
%%% apply threshold to find patches
superthresh = imroiblurred>=thresh;
subthresh = imroiblurred<thresh;
[bnds_allsizes, labelmat_allsizes, nblobsraw] = bwboundaries(superthresh,'noholes');%% find patches
patchdata.thresh = thresh;

%%% optionally include non-roi pixels that are completely contained within patches as part of the surrounding patch
if pars.include_interior_nonroi_in_patches
    [~, lm, ~, adj] = bwboundaries(logical(labelmat_allsizes));
    [inner_nan, surrounding_patch] = find(adj);
    for i = 1:length(inner_nan)
        superthresh(lm==inner_nan(i)) = true; %% fill in nan holes
        [bnds_allsizes, labelmat_allsizes, nblobsraw] = bwboundaries(superthresh,'noholes');%% find patches again with nan holes filled
        if strcmp(pars.threshMode,'intensityQuantiles') % fill in nan holes in quantiles (highest-level quantile touching the nans)
            quantfound = false;
            iquant = nQuantiles;
            while ~quantfound
                if any(find(quantile_table.edges{iquant} & lm==inner_nan(i))) % if this nan hole touches this quantile
                    quantile_table.quantimage{iquant}(lm==inner_nan(i)) = true; % add nans to this quantile
                    quantile_table.edges{iquant} = edge(quantile_table.quantimage{iquant}); % get edges of this quantile again
                    quantfound = true;
                elseif iquant > 1
                    iquant = iquant-1;
                else
                    quantfound = true; 
                end
            end
        end
    end
end

% get image of all edges, for illustration only; bottom level = outside of Roi
quant_edges_img = zeros(size(im)); % do not use for quantification, illustration only; bottom level = outside of Roi
for iquant = 1:nQuantiles
    quant_edges_img(quantile_table.edges{iquant}) = iquant;
end
patchdata.quant_edges_img = quant_edges_img; % do not use for quantification, illustration only; bottom level = outside of Roi


    




