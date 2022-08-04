%%% find cells from black-white image, count cells in each quantile from an 
%%% analyzed patch image, perform chi square test to look for signficant
%%% different in cell density per quantile
% [chisq_pval, quant_cells_table] = chi_square_quantiles_cells(patchdata, cellsimage)
%%% patchdata == output from findpatches.m; must include quantile analysis
%%% newroi (optional) == bw mask for specifying which region of the image should be analyzed
%%%     -quant images from patchdata must include all newroi pixels
%%%%%% last updated 18/9/27 on thermaltake

function [chisq_pval, quant_cells_table, cell_centers_img] = chi_square_quantiles_cells(patchdata, cellsimage, newroi)

cellsimage = loadbw(cellsimage);
[cellCentersInds_raw, cellCentersInds_roiRelative, cell_centers_img] = find_shape_centers(cellsimage, patchdata.roi);
nQuantiles = patchdata.pars.nQuantiles;
if exist('newroi','var') % if we are going to cut out a section of the roi that was originally used to find quantiles/patches
    newroi = loadbw(newroi);
    for iquant = 1:nQuantiles
        patchdata.quantile_table.quantimage{iquant} = patchdata.quantile_table.quantimage{iquant} & newroi; % apply newroi mask
        patchdata.quantile_table.pixyx{iquant} = []; % get quant pix coordinates again
        [patchdata.quantile_table.pixyx{iquant}(:,1), patchdata.quantile_table.pixyx{iquant}(:,2)] = find(patchdata.quantile_table.quantimage{iquant});
    end
end
quant_cells_table = patchdata.quantile_table(:,'quantimage');
quant_cells_table.ncells = NaN(nQuantiles,1);
quant_cells_table.cells_per_sqmm = NaN(nQuantiles,1);
xyMicronsPerPixel = umPerPix(patchdata.zoom,patchdata.scope);
sqMMPerPixel = xyMicronsPerPixel^2 / 1e6; % area of square pixel in sq mm

for iquant = 1:nQuantiles
    npix = size(patchdata.quantile_table.pixyx{iquant},1); % number of pix in this quantile
    nSqMM = sqMMPerPixel * npix; % area of this quantile in square mm
    quant_cells_table.ncells(iquant) = nnz(quant_cells_table.quantimage{iquant} & cell_centers_img); % get number of cells falling in this quantile
    quant_cells_table.cells_per_sqmm(iquant) = quant_cells_table.ncells(iquant) / nSqMM; 
end

% chi square test
freqs = quant_cells_table.ncells;
noutcomes = length(freqs);
degfr = noutcomes-1;
expected = sum(freqs) / noutcomes; % frequencies if they were evenly distributed across outcomes
devstat = (freqs-expected).^2 / expected; % calculate deviations from expected for each outcome
chistat = sum(devstat); % chi2 test statistic
chisq_pval = chi2cdf(chistat,degfr,'upper'); % get pval

