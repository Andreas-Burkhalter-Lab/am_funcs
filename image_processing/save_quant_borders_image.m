%%%% SAVE_QUANT_BORDERS_IMAGE
%%%% save image of quantiles from a patchdata struct, with alpha masking, that can be overlaid with photoshop images
%
% save_quant_borders_image(patchdata,ops)
%   ops struct includes fields: quants, border_thickness_pix, savepars
%   savepars is the options structure for savetif.m; important options include: conv_factor, zero_transparent
%
% updated 20/10/15 on thermaltake

function save_quant_borders_image(patchdata,ops)

%% params
ops = vardefault('ops',struct);
ops.savepars = field_default(ops,'savepars',struct);
ops.quants = field_default(ops,'quants',[ 5]); % if threshMode==intensityQuantiles, draw the borders between these quants and the quants directly below them
ops.border_thickness_pix = field_default(ops,'border_thickness_pix',5); % contour thickness in pixels
ops.colored_borders = field_default(ops,'colored_borders',false); % if false, image will be logical; if true, image will be RGBA... quantiles only
    ops.border_colors = field_default(ops,'border_colors',{[1 0 0], [0 1 0], [0 0 1]}); % colors to assign to each quantile border
    ops.color_intensity_multiplier = field_default(ops,'color_intensity_multiplier',2^16); % for RGB image, multiply intensity by this value
ops.savepars.conv_factor = field_default(ops.savepars,'conv_factor',16^3.5); % brightness
ops.savepars.zero_transparent = field_default(ops.savepars,'zero_transparent',1); % non-border pixels will be transparent

%%
if isfield(patchdata,'roifile') && ~isempty(patchdata.roifile) && ischar(patchdata.roifile)
    save_prepend = strrep(getfname(patchdata.roifile), 'roi', '');
else
    save_prepend = '_'; 
end

if strcmp(patchdata.pars.threshMode, 'intensityQuantiles') % draw quantiles if this patchfinding method was used
    edges_img = patchdata.quant_edges_img;
    keeppix = false(size(edges_img));
    for iquant = 1:length(ops.quants)
        quantval = ops.quants(iquant);
        keeppix = keeppix | edges_img==quantval;
    end
    edges_img(~keeppix) = 0; % clear edges values not listed in quants
    savename = [save_prepend, 'quantborders' num2str(ops.quants)];
    

elseif ~strcmp(patchdata.pars.threshMode, 'intensityQuantiles') % if quantiles weren't found, just draw patch borders
    edges_img = edge(patchdata.patchimage);
    savename = [save_prepend, 'patchborders'];
end

%%% thicken pixels by added shifted version of the image; alternate directions to shift
edges_img_thick = edges_img;
n_lines_to_add = ops.border_thickness_pix-1; 
for iline = 1:n_lines_to_add
    shift_n_pix = ceil(iline/4);
    direction_to_shift = mod(iline,4);
    switch direction_to_shift
        case 1 % shift up
            shifted_img = circshift(edges_img,[-shift_n_pix, 0]);
            shifted_img(end,:) = 0; % don't overflow to other side of image
        case 2 % shift right 
            shifted_img = circshift(edges_img,[0, shift_n_pix]);
            shifted_img(:,1) = 0; % don't overflow to other side of image
        case 3 % shift down
            shifted_img = circshift(edges_img,[shift_n_pix, 0]);
            shifted_img(1,:) = 0; % don't overflow to other side of image
        case 0 % shift left
            shifted_img = circshift(edges_img,[0, shift_n_pix]);
            shifted_img(:,end) = 0; % don't overflow to other side of image
    end
    pixels_to_modify = logical(shifted_img);
    edges_img_thick(pixels_to_modify) = shifted_img(find(shifted_img)); % use replacement rather than summation because we don't want to create new pixel values
end

% if borders need to be colored
if ops.colored_borders &&  strcmp(patchdata.pars.threshMode, 'intensityQuantiles')
    edges_img_thick_1d = edges_img_thick; 
    edges_img_thick = zeros(size(edges_img_thick,1), size(edges_img_thick,2), 3); % 3D for RGB; alpha layer will be added in savetif.m
    for iquant = 1:length(ops.quants)
        quantval = ops.quants(iquant);
        this_quant_color = ops.border_colors{iquant}; % RGB color to assign
        this_quant_border_img = double(edges_img_thick_1d==quantval); % pixels to color
        for dim = 1:3
            edges_img_thick(:,:,dim) = edges_img_thick(:,:,dim) + this_quant_color(dim)*this_quant_border_img; % assign color to this RGB dimension
        end
    end
%     edges_img_thick(:,:,4) = logical(edges_img_thick_1d); % add alpha; all zero pixels are transparent
    edges_img_thick = ops.color_intensity_multiplier * edges_img_thick; % set intensity of output image
end

savetif(edges_img_thick, savename, ops.savepars)
fprintf(['Saved borders image ' savename '.tif\n'])