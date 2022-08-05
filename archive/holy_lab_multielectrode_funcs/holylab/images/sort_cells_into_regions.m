function cells_in_regions = sort_cells_into_regions(cell_points,regions)
% SORT_CELLS_INTO_REGIONS tests whether each cell_point lies within the
%                               polygonal regions described for each image
%                               within the regions array using the inpolygon fxn.
%                               It then sums the cells within each region,
%                               and returns cells_in_regions, which is a
%                               cell array the same size as regions
%                               
% See also CELLCOUNT_FLAT, INPOLYGON
%
 
% Copyright 2007 Julian P. Meeks

% determine size of input "regions"
n_images = size(regions, 2);
for idx = 1:n_images
    n_regions(idx) = size(regions{idx},2);
    n_cellpoints(idx) = size(cell_points{idx},2);
end

% change coordinate points stored in regions into arrays of x and y pts
region_x_pts = cell(size(regions));
region_y_pts = cell(size(regions));
cellpoints_x_pts = cell(size(regions));
cellpoints_y_pts = cell(size(regions));

for idx_img = 1:n_images
    for idx_reg = 1:n_regions(idx_img)
        for idx_point = 1:length(regions{idx_img}{idx_reg})
            region_x_pts{idx_img}{idx_reg}(idx_point)=regions{idx_img}{idx_reg}{idx_point}(1);
            region_y_pts{idx_img}{idx_reg}(idx_point)=regions{idx_img}{idx_reg}{idx_point}(2);
        end
    end
end
for idx_img = 1:n_images
    for idx_cpts = 1:n_cellpoints(idx_img)
        cellpoints_x_pts{idx_img}(idx_cpts) = cell_points{idx_img}{idx_cpts}(1);
        cellpoints_y_pts{idx_img}(idx_cpts) = cell_points{idx_img}{idx_cpts}(2);
    end
end

% calculate areas
sorted_cells = cell(size(cell_points));
cells_in_regions = cell(size(regions));

for idx_img = 1:n_images
    for idx_reg = 1:n_regions(idx_img)
        for idx_cpts = 1:n_cellpoints(idx_img)
            sorted_cells{idx_img}(idx_cpts)= ...
                    inpolygon(cellpoints_x_pts{idx_img}(idx_cpts),...
                              cellpoints_y_pts{idx_img}(idx_cpts),...
                              region_x_pts{idx_img}{idx_reg},...
                              region_y_pts{idx_img}{idx_reg});
        end
        cells_in_regions{idx_img}(idx_reg) = sum(sorted_cells{idx_img});
    end
end
    
end