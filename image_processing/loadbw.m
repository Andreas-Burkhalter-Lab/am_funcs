%%%% load bw image indicating region of interest
%
% [img imfile] = loadbw(imfileIn, ops)
% imfileIn can be either filename or image
% ops structure includes:
%   logical_zero_value = value(s) in input to be used as 'false' value; all other values are 'true' (default 0)
%   allow_multiple_nonzero: if false, generate error if there are multiple nonzero values (default false)
%   reverse_black_white: if true, reverse true and false (after applying logical_zero_value) 
%%%%% updated 2020/8/14 on thermaltake

function [img imfile] = loadbw(imfileIn, ops)

ops = vardefault('ops',struct);
ops.logical_zero_value = field_default(ops,'logical_zero_value',0); % if not specified, assume that 'false' value is zero
ops.allow_multiple_nonzero = field_default(ops,'allow_multiple_nonzero',false); % default to requiring only 1 nonzero value
ops.reverse_black_white = field_default(ops,'reverse_black_white',false); % default to not inverting 
if ischar(imfileIn)
    input_is_filename = 1; 
    imfile = imfileIn;
    img = imread(imfileIn);
else
    input_is_filename = 0; 
    imfile = '';
    img = imfileIn;
end
    
img = img(:,:,1); % in case image is rgb, take only first channel

% convert the false value override to zero
if any(ops.logical_zero_value ~= 0)
    old_zeros = img == 0;
    new_zeros = false(size(img)); 
    for i = 1:length(ops.logical_zero_value) % check for elements that match any zero values
        new_zeros = new_zeros | [img == ops.logical_zero_value(i)];
    end
    img(old_zeros) = ops.logical_zero_value(1); % old zeros become 'true'
    img(new_zeros) = 0; % zero override values become 'false'
end

% optionally invert black and white
if ops.reverse_black_white
    img = ~img;
end

if ~ops.allow_multiple_nonzero && numel(unique(img(img~=0))) ~= 1 % check for multiple nonzero vals
    error('More or less than 1 unique non-zero pixel value found in ROI file.')
end
img = img>0;

if input_is_filename %%% if image is from file, check for border artifacts
    if any(find(img(:,1))) || any(find(img(:,end))) || any(find(img(1,:))) || any(find(img(end,:)))
        f = figure('units','normalized','outerposition',[0 0 1 1]);
        imagesc(img)
        pause(.1)
        commandwindow
        go_on = input(['Warning: BW image ' imfile ' has nonzero pixels along a border. Enter ''y'' to continue.'],'s');
        close(f)
        if ~strcmp(go_on,'y')
            error('quitting loadbw')
        end
    end
end