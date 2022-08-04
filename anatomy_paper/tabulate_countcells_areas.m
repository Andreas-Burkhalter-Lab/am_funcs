%%% get cell counts and densities in rois
% input filetable must have variables 'cellsfile','roifile','zoom','scope'
% updated 2/17/18

function [ftable_out, lamtable_totalcells, lamtable_cellsPerSqmm] = countcells(filetable)

nfiles = height(filetable);
ftable_out = filetable;
ftable_out.cellsInRoiYX = cell(nfiles,1);
ftable_out.ncells = NaN(nfiles,1);
ftable_out.roiAreaSqMM = NaN(nfiles,1);
ftable_out.cellsPerSqMM = NaN(nfiles,1);

% get cell counts and density for each file
for i = 1:nfiles
    clear cellsInRoiYX cellsInRoiInds cellcentersInds xyPixelsPerUM sqMMPerPixel
    xyPixelsPerUM = pixPerUm(filetable.zoom(i),filetable.scope{i}); 
    sqMMPerPixel = xyPixelsPerUM^-2 / 10e6; % area of square pixel in sq mm
    
    
    roi = loadbw(ftable_out.roifile{i});
    nRoiPixels = length(find(roi));
    ftable_out.roiAreaSqMM(i) = nRoiPixels * sqMMPerPixel;
    ftable_out.cellsfile{i}
    cellimage = loadbw(ftable_out.cellsfile{i});
    [bounds, labelmat, ncells_total, adjc] = bwboundaries(cellimage);
    cellcentersYX = NaN(ncells_total,2);
    for j = 1:ncells_total
        [y x] = find(labelmat == j);
        cellcentersYX(j,1) = round(mean(y));
        cellcentersYX(j,2) = round(mean(x));
    end
    roiInds = find(roi);
    cellcentersInds = sub2ind(size(roi),cellcentersYX(:,1),cellcentersYX(:,2));
    cellsInRoiInds = intersect(cellcentersInds,roiInds);
    [cellsInRoiYX(:,1) cellsInRoiYX(:,2)] = ind2sub(size(roi),cellsInRoiInds);
    ftable_out.ncells(i) = length(cellsInRoiInds);
    ftable_out.cellsPerSqMM(i) = ftable_out.ncells(i) / ftable_out.roiAreaSqMM(i);
end



% tabulate cell counts and density by area and layer
if exist('lamtable_totalcells','var')
    nn = NaN(length(unique(cellcount.sec)),1);
    lamtable_totalcells = table(nn,nn,nn,nn,nn,nn,'VariableNames',{'sec','pora','por','p','lm','li'});
    lamtable_totalcells.sec = unique(cellcount.sec);
    lamtable_cellsPerSqmm = lamtable_totalcells;

    q = find(strcmp(cellcount.area,'pora'));
    lamtable_totalcells.pora = cellcount.ncells(q);
    lamtable_cellsPerSqmm.pora = cellcount.cellsPerSqMM(q);
    q = find(strcmp(cellcount.area,'por'));
    lamtable_totalcells.por = cellcount.ncells(q);
    lamtable_cellsPerSqmm.por = cellcount.cellsPerSqMM(q);
    q = find(strcmp(cellcount.area,'p'));
    lamtable_totalcells.p = cellcount.ncells(q);
    lamtable_cellsPerSqmm.p = cellcount.cellsPerSqMM(q);
    q = find(strcmp(cellcount.area,'lm'));
    lamtable_totalcells.lm = cellcount.ncells(q);
    lamtable_cellsPerSqmm.lm = cellcount.cellsPerSqMM(q);
    q = find(strcmp(cellcount.area,'li'));
    lamtable_totalcells.li = cellcount.ncells(q);
    lamtable_cellsPerSqmm.li = cellcount.cellsPerSqMM(q);
end