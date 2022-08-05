function region_areas = calculate_region_area(regions)
% CALCULATE_REGION_AREA calculates the area of a given region from
%                               hist_analysis field in .xdb database
%                               and returns region_areas, which is a
%                               cell array of the same size as regions
%                               containing one integer value per region
%
% See also CELLCOUNT_FLAT, POLYAREA
%
 
% Copyright 2007 Julian P. Meeks

% determine size of input "regions"
n_images = size(regions, 2);
for idx = 1:n_images
    n_regions(idx) = size(regions{idx},2);
end

% change coordinate points stored in regions into arrays of x and y pts
x_pts = cell(size(regions));
y_pts = cell(size(regions));

for idx_img = 1:n_images
    for idx_reg = 1:n_regions(n_images)
        for idx_point = 1:length(regions{idx_img}{idx_reg})
            x_pts{idx_img}{idx_reg}(idx_point)=regions{idx_img}{idx_reg}{idx_point}(1);
            y_pts{idx_img}{idx_reg}(idx_point)=regions{idx_img}{idx_reg}{idx_point}(2);
        end
    end
end

% calculate areas
region_areas = cell(size(regions));

for idx_img = 1:n_images
    for idx_reg = 1:n_regions(n_images)
        region_areas{idx_img}{idx_reg}= double(...
                    polyarea(x_pts{idx_img}{idx_reg},y_pts{idx_img}{idx_reg}));
    end
end

end
