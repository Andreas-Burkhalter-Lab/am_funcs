 %%% find max number of patches
 
 function [ analysis_data patchtable] = find_most_patches( imfile, roifile, zoom, scope, pars )

plotting = 0;
output_table_images = 0;  

if ~exist('pars','var') || ~isfield(pars,'threshMode') 
    pars.threshMode = 'fractionOfMax';
end
if ~exist('pars','var') || ~isfield(pars,'threshFracs') 
    threshFracs = 0.8:0.002:0.96;
else
    threshFracs = pars.threshFracs;
end
if ~exist('scope','var')
    scope = [];
end

patchtable = table(threshFracs','VariableNames',{'threshFrac'});

for i = 1:length(threshFracs)
    pars.threshFraction = threshFracs(i);
    p = findpatches(imfile,roifile,zoom,scope,pars);
    patchtable.nlarge(i) = length(p.bnds_large);
    patchtable.nsmall(i) = length(p.bnds_small);
    patchtable.bnds_large{i} = p.bnds_large;
    patchtable.bnds_large_image{i} = p.bnds_large_image;
    patchtable.bnds_small_image{i} = p.bnds_small_image;
    patchtable.labelmat_large{i} = p.labelmat_large; 
end
[analysis_data.maxpatches ind] = max(patchtable.nlarge);
analysis_data.maxPatchesThreshFrac = patchtable.threshFrac(ind);

analysis_data.roi = p.roi;
npixroi = numel(find(analysis_data.roi)); 
analysis_data.roi_area_sqmm = npixroi/[pixPerUm(zoom,scope)]^2 * 0.001^2;
analysis_data.patches_per_sqmm = analysis_data.maxpatches / analysis_data.roi_area_sqmm;
analysis_data.threshMode = p.threshMode;
analysis_data.minAreaPatch_squm = p.minAreaPatch_squm; % min number of pix a blob must contain to be considered a patch
analysis_data.diskblurradius_um = p.diskblurradius_um; 
analysis_data.imblurred = p.imblurred;
analysis_data.imroiblurred = p.imroiblurred;
analysis_data.labelmat_large = patchtable.labelmat_large{ind};
analysis_data.bnds = patchtable.bnds_large{ind};
analysis_data.bnds_large_image = patchtable.bnds_large_image{ind};
analysis_data.bnds_small_image = patchtable.bnds_small_image{ind};
analysis_data.patchimage = logical(analysis_data.labelmat_large);
analysis_data.pixPerUm = pixPerUm(zoom,scope);
analysis_data.zoom = zoom;
analysis_data.scope = scope;


if ~output_table_images
    patchtable.bnds_large = [];
    patchtable.bnds_large_image = [];
    patchtable.bnds_small_image = [];
    patchtable.labelmat_large = [];
end    

if plotting
    figure('units','normalized','outerposition',[0 0 1 1])
    plot(threshFracs,patchtable.nlarge)
    xlabel('threshfrac')
    ylabel('npatches')
    title(imfile)
    commandwindow
end