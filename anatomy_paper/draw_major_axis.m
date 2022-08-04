%%% draw_major_axis: draw major axis of area defined by white shapes in BW image
% intended for use with only one white object in input image
%
% major_axis = draw_major_axis(imfile)
%
% updated 2019/12/30

function major_axis_image = draw_major_axis(imfile)

 save_output = 1; % save major axis image as tif file
    append_tag = '_axis'; % append this to the input image filename
% % %  line_thickness_pix = 3; 
 
im = loadbw(imfile);
savename = [getfname(imfile), append_tag];

% find major axis
area_centroid = regionprops(im,'Centroid');
    centroid_x = area_centroid.Centroid(1);
    centroid_y = area_centroid.Centroid(2);
orientation = regionprops(im,'Orientation'); % major axis angles
    orientation = orientation.Orientation;
axis_slope = -tan(deg2rad(orientation)); % slope of major axis line to draw
imagelength = size(im,2);
imageheight = size(im,1);
axis_x = 1:imagelength; % x coordinates of axis to draw
y_intercept = centroid_y - axis_slope*centroid_x;
axis_y = axis_slope*axis_x + y_intercept; % y coordinates of axis to draw
axis_x = round( axis_x(axis_y < imageheight & axis_y > 0) ); % truncate x axis coordinates with image bounds
axis_y = round( axis_y(axis_y < imageheight & axis_y > 0) ); % truncate y axis coordinates with image bounds

% draw major axis
major_axis_image = false(size(im));
inds = sub2ind(size(major_axis_image), axis_y, axis_x);
major_axis_image(inds) = true;
% imagesc(major_axis_image)
savetif(major_axis_image, savename)