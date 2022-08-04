%% this version does not automatically create a pixel buffer if the roi includes border pixels
%%

%%%% find threshold intensity for patches
%%% updated 19/12/13 on thermaltake

%%% blurring
H = fspecial('disk',diskblurradius_pix); %%% create filter; arguments set blur radius
imroi = im;
imroi = double(imroi);  % must be double to have nans
imroi(~roi) = NaN;
patchdata.imroi  = imroi; 
if pars.blur_using_only_roi
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
        patchimage = false(size(im)); % patch pixels will be changed to true
        % divide image into intensity quantiles
        for iquant = 1:nQuantiles
            if iquant == 1
                quantile_table.quantimage{iquant} = imroiblurred <= quantvals(iquant) & roi;
            else % eliminate overlaps between quantiles
                quantile_table.quantimage{iquant} = imroiblurred <= quantvals(iquant) & imroiblurred > quantvals(iquant-1);
                quantile_table.pixIntensThresh(iquant) = quantvals(iquant-1); %% minimum pix intens to equal or exceed to be included in this quant group
            end
            [quantile_table.pixyx{iquant}(:,1), quantile_table.pixyx{iquant}(:,2)] = find(quantile_table.quantimage{iquant});
            quant_levels_img = quant_levels_img + iquant*double(quantile_table.quantimage{iquant});
            quantile_table.edges{iquant} = edge(quantile_table.quantimage{iquant}); % get edges of this quantile
            if any(iquant == pars.patchQuantiles) % if this quantile is to be labeled as patches
                patchimage = patchimage | quantile_table.quantimage{iquant}; % add this quant to patches image
            end
        end
%         quant_levels_img(~roi) = NaN; % make sure that out-of-roi pixels aren't included in a quantile yet
        patchdata.quantile_table = quantile_table;
        patchdata.quant_levels_img = quant_levels_img;  
    case 'fractionOfPopulation'
        roisorted = sort(imroivals);
        threshpoint = max([1, length(roisorted) - round(pars.fracPixSuperthresh * length(roisorted))]);
        thresh = roisorted(threshpoint); 
        patchdata.threshpoint = threshpoint;
        patchdata.threshFraction_fractionOfMax = double(thresh)/double(max(roisorted)); % convert for comparison with other threshmode
        threshinc = maxInRoi * pars.raisethresh_increment;
    case 'fractionOfMaxMinusMin'
        thresh =  minInRoi + [maxInRoi-minInRoi] * pars.threshFraction;
        threshinc = [maxInRoi-minInRoi] * pars.raisethresh_increment;
    case 'fractionOfMax'
        thresh =  maxInRoi * threshFraction;
        threshinc = maxInRoi * pars.raisethresh_increment;
end
    
%%% apply threshold and find patches
if ~strcmp(pars.threshMode,'intensityQuantiles') %% if a single threshold is to be found
    patchimage = imroiblurred>=thresh;
    patchdata.thresh = thresh;
    [bnds_allsizes, labelmat_allsizes, nblobsraw, ~] = bwboundaries(patchimage,'noholes');%% find patches
end


%%% optionally include non-roi pixels that are completely contained within a quantile/patch as part of the containing region
%%%%%%%  ... allows us to analyze cells that fall on blood vessels
if pars.include_interior_nonroi_in_roi && nnz(~roi)>0 % only run this step if there are pixels marked as non-roi
    if any(find(roi(:,1))) || any(find(roi(:,end))) || any(find(roi(1,:))) || any(find(roi(end,:))) % if roi pixels are on the border
        checkinput = input('\npars.include_interior_nonroi_in_roi==true is only designed to work with a fully non-roi border; continue? ','s');
        if ~strcmp(checkinput,'y')
            error('quitting function')
        end
        warning('check patch border assignment for errors')
    end
    if strcmp(pars.threshMode,'intensityQuantiles') % if we need to identify a quantile for the nonroi spot
        [nan_bnds, nan_lm] = bwboundaries(roi);
        for iblob = 1:max(max(nan_lm))
            this_blob_im = nan_lm==iblob;
            this_blob_im = this_blob_im & quant_levels_img==0; % only look at the part of the non-roi blob that has not yet been assigned to a quantile
            if any(any(this_blob_im & ~roi)) % if this blob was marked as outside the roi
                clear blob_yx
                [blob_yx(:,1), blob_yx(:,2)] = find(this_blob_im);
                blob_with_outer_bnds = [blob_yx + [1 1]; blob_yx + [-1 1]; blob_yx + [1 -1]; blob_yx + [-1 -1]];
                blob_with_outer_bnds_inds = sub2ind(size(im),blob_with_outer_bnds(:,1),blob_with_outer_bnds(:,2));
                bording_quant_vals = quant_levels_img(blob_with_outer_bnds_inds);
                mode_bordering_quantile = mode(bording_quant_vals(bording_quant_vals>0)); % finding most frequent quantile touching this non-roi blob
                quant_levels_img(this_blob_im) = mode_bordering_quantile;
                quantile_table.quantimage{mode_bordering_quantile}(this_blob_im) = true; % add former nans to this quantile
                quantile_table.edges{mode_bordering_quantile} = edge(quantile_table.quantimage{mode_bordering_quantile}); % redraw edges of the quantile the blob is joining
                if any(mode_bordering_quantile == pars.patchQuantiles) % if this quantile is to be labeled as patches
                    patchimage = patchimage | quantile_table.quantimage{mode_bordering_quantile}; % redraw patches
                end
            end
        end
        for iquant = 1:nQuantiles  % update pix coordinates
            quantile_table.pixyx{iquant} = []; % get quant pix coordinates again
            [quantile_table.pixyx{iquant}(:,1), quantile_table.pixyx{iquant}(:,2)] = find(quantile_table.quantimage{iquant});
        end
        patchdata.quantile_table = quantile_table;
        patchdata.quant_levels_img = quant_levels_img;  
        
    %%% optionally include non-roi pixels that are completely contained within a quantile/patch as part of the containing region
    elseif ~strcmp(pars.threshMode,'intensityQuantiles') % if we don't need to identify a quantile for the nonroi spot
%         error('check that this code section still makes sense')
        [~, lm, ~, adj] = bwboundaries(logical(labelmat_allsizes));
        [inner_blob, surrounding_patch] = find(adj);
        for i = 1:length(inner_blob)
            patchimage(lm==inner_blob(i)) = true; %% fill in hole
            [bnds_allsizes, labelmat_allsizes, nblobsraw] = bwboundaries(patchimage,'noholes');%% find patches again with nan holes filled
        end
    end
end

% get image of all edges, for illustration only; bottom level = outside of Roi
if strcmp(pars.threshMode,'intensityQuantiles')
    quant_edges_img = zeros(size(im)); % do not use for quantification, illustration only; bottom level = outside of Roi
    for iquant = 1:nQuantiles
        quant_edges_img(quantile_table.edges{iquant}) = iquant;
    end
    patchdata.quant_edges_img = quant_edges_img; % do not use for quantification, illustration only; bottom level = outside of Roi
end

    




