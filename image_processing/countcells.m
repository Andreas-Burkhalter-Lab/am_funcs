%%% get cell counts and densities in rois
% [ftable_out] = countcells(filetable, table_rows_, roi)
%
%%%%%%%%% input filetable (table or spreadsheet filename) must have variables 'cellsfile','zoom','scope'
% optionally include vector 'excel_rows' to specify which rows to analyze (assumes top filler row)
%%%%%%%%% either include 'roi' input to apply to all cell images or include
%%%%%%%%% variable 'roifile' in filetable to apply individually to each
%%%%%%%%% cell image
%   optionally include variable 'regfile' in filetable for warping images 
%   (see assign_fiducials.m)
%
% updated 2020/5/10 on thermaltake

function [ftable_out] = countcells(filetable, excel_rows, roi)

if ischar(filetable) % open excel file containing filetable
    ftable_out = readtable(filetable);
elseif istable(filetable) % if input is table variable
    ftable_out = filetable;
end
excel_rows = vardefault('excel_rows',2:1+height(ftable_out)); % default to using all table rows
ftable_out = ftable_out(excel_rows-1,:); % only include the specified case rows; assumes top excel row was discarded by readtable
nfiles = height(ftable_out);

% % % ftable_out.ncells = NaN(nfiles,1);
ftable_out.cellsPerSqMM = NaN(nfiles,1);
ftable_out.cellsYX = cell(nfiles,1);
if any(strcmp('regfile',ftable_out.Properties.VariableNames)) % if the cell images need to be warped
    ftable_out.cellimage_prereg = cell(nfiles,1);
    fprintf('Performing cell image warping.')
end
 

% get cell counts and density for each file
for i = 1:nfiles
    %%% if the cell was already assigned a cellcount in the excel file (~isnan), don't analyze it here
    if isnan(ftable_out.ncells(i)) 
        clear cellsYX cellsInRoiInds cellcentersInds xyPixelsPerUM sqMMPerPixel
        if ~any(strcmp('scope',ftable_out.Properties.VariableNames))
            scope = 'epimicro';
        else
            scope = ftable_out.scope{i};
        end
        xyPixelsPerUM = pixPerUm(ftable_out.zoom(i),scope); 
        sqMMPerPixel = xyPixelsPerUM^-2 / 10e6; % area of square pixel in sq mm

        % get roi in which to count cells
        if exist('roi','var')
            if any(strcmp('roifile',ftable_out.Properties.VariableNames))
                error('roi specified twice')
            end
            [roi_thisfile roifile] = loadbw(roi);
            nRoiPixels = length(find(roi)); 
            ftable_out.Properties.UserData.roiAreaSqMM = nRoiPixels * sqMMPerPixel;
            roiAreaSqMM_thisfile = ftable_out.Properties.UserData.roiAreaSqMM;
        else    
            roi_thisfile = loadbw(ftable_out.roifile{i});   
            nRoiPixels = length(find(roi_thisfile)); 
            ftable_out.roiAreaSqMM(i) = nRoiPixels * sqMMPerPixel;
            roiAreaSqMM_thisfile = ftable_out.roiAreaSqMM(i);
        end

        ftable_out.cellsfile{i}; %% take off semicolon to track progress w/ command line output
        cellimage = loadbw(ftable_out.cellsfile{i});
        %%%%% first get cell centers from prereg image, then warp image of cell
        %%%%% centers, to reduce chances of creating cell overlaps through
        %%%%% warping
        if any(strcmp('regfile',ftable_out.Properties.VariableNames)) % if the cell images need to be warped
            cellimage_prereg = cellimage;   
            ftable_out.cellimage_prereg{i} = cellimage_prereg;
            reg_struct = load(ftable_out.regfile{i}); % get params for warping
            if any(size(cellimage_prereg) ~= size(reg_struct.movingimage))
                error('moving image dimension mismatch')
            end
            [~, labelmat_prereg, ncells_total_prereg, ~] = bwboundaries(cellimage); % get cell locs from prereg image
            cellcenters_im_prereg = false(size(labelmat_prereg));
            % get cells from prereg image
            for icell = 1:ncells_total_prereg
                [y x] = find(labelmat_prereg == icell);
                centery = round(mean(y));
                centerx = round(mean(x));
                cellcenters_im_prereg(y,x) = true; % label cell center
            end
            cellimage = imwarp(cellcenters_im_prereg,reg_struct.geotransform,'OutputView',reg_struct.worldref);
        end

        % get cell centers
        [bounds, labelmat, ncells_total, adjc] = bwboundaries(cellimage);
        cellcentersYX = NaN(ncells_total,2);
        for icell = 1:ncells_total
            [y x] = find(labelmat == icell);
            cellcentersYX(icell,1) = round(mean(y));
            cellcentersYX(icell,2) = round(mean(x));
        end
        roiInds = find(roi_thisfile);
        cellcentersInds = sub2ind(size(roi_thisfile),cellcentersYX(:,1),cellcentersYX(:,2));
        cellsInRoiInds = intersect(cellcentersInds,roiInds);
        [cellsYX(:,1) cellsYX(:,2)] = ind2sub(size(roi_thisfile),cellsInRoiInds);
        ftable_out.cellimage{i} = cellimage; 
        ftable_out.ncells(i) = length(cellsInRoiInds);
        ftable_out.cellsYX{i} = cellsYX;
        ftable_out.cellsPerSqMM(i) = ftable_out.ncells(i) / roiAreaSqMM_thisfile;
    end
end

% rearrange variables
ftable_out = movevars(ftable_out, 'cellsPerSqMM', 'Before', 'scope');
ftable_out = movevars(ftable_out, 'ncells', 'Before', 'cellsPerSqMM');