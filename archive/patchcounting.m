 % input has variables: dir, m2file, roifile, zoom
 
 function [patchcount] = patchcounting(filetable)
 
threshFracs = 0.8:0.002:0.96; % see find_most_patches.m
 
nfiles = height(filetable);
patchcount = filetable;
patchcount.mostpatches_table = cell(nfiles,1);
 
 for i = [1:nfiles]
    [p a] = find_most_patches([patchcount.dir{i} filesep patchcount.m2file{i}],[ patchcount.dir{i} filesep patchcount.roifile{i}], threshFracs);
    patchcount.mostpatches_table{i} = p;
    patchcount.mostpatches_analysis{i} = a;
    patchcount.npatches(i) = patchcount.mostpatches_analysis{i}.maxpatches;
    npixroi = numel(find(patchcount.mostpatches_analysis{i}.roi)); 
    patchcount.por_area_sqmm(i) = npixroi/[zoomConvert(patchcount.zoom(i))]^2 * 0.001^2;
    patchcount.patches_per_sqmm(i) = patchcount.npatches(i) / patchcount.por_area_sqmm(i);
 end