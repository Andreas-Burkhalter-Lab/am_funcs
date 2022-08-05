function mask = stack_automask(img_in, ops)
% stack_automask takes an input img and creates a custom 3d mask
% Syntax: mask = stack_automask(img_in, ops)
% Inputs:
%    img_in is a grayscale image stack (2d or 3d) 
%    ops is a structure containing the following optional fields
%           .deadframes: a one or two-vector containing the number of masked
%                        frames at the beginning and end of the stack
%                        (default: 4)
%           .deadperim:  a one, two, or four-vector containing the number of masked
%                        pixels at the x and y edges of the stack (default: 15)
%                        i.e. [4] masks out 4 pixels on all edges
%                             [4 5] masks out 4 x-dim pixels and 5 y-dim
%                             [4 10 5 11] uses values to mask [left,right,bottom,top]
%           .floor:        an integer minimum pixel threshold value 
%                          (default []) if empty, no floor will be applied
%           .nbins:      : an integer value containing the number of bins
%                          to use for image histogram (default 128)
%           .filter:     : an image filter (see fspecial) to smooth mask
%                          edges (default: fspecial('gaussian',16,4))
%           .scalefactor : a scalar multiplier of the median point (default 1)
%                          values in the range 0 < scalefactor < 1 reduce
%                             the threshold towards the floor_offset value
%                          values > 1 raise the intensity towards the
%                             maximal value (clipped at max 14-bit intensity for now)
% Outputs:
%    mask is a binary matrix of the same size as img_in with values of 1 at
%         kept pixels and 0 at discarded pixels
%
% See also register_multigrid_options, fspecial

% Copyright 2010 Julian P. Meeks

%% check inputs
stack_size = size(img_in);
stack_size(end+1:3) = 1;

if ~exist('ops','var');
    ops = struct;
elseif ~isstruct(ops)
    warning('''ops'' was not supplied as a structure variable as intended.  using default values instead');
    ops = struct;
end
ops = default(ops, 'deadframes',4);
ops = default(ops, 'deadperim',15);
ops = default(ops, 'floor',[]);
ops = default(ops, 'nbins',128);
ops = default(ops, 'filter',fspecial('gaussian',16,4));
ops = default(ops, 'scalefactor',1);

if length(ops.deadframes) < 2
    ops.deadframes = repmat(ops.deadframes,1,2);
end
if length(ops.deadperim) == 1
    ops.deadperim = repmat(ops.deadperim,1,4);
elseif length(ops.deadperim) == 2
    temp(1:2) = ops.deadperim(1);
    temp(3:4) = ops.deadperim(2);
    ops.deadperim = temp;
end

%% set "floor" of intensities to be masked
ibins = 0:(2^14-1)/(ops.nbins-1):2^14-1;
mask = zeros(stack_size(1:3));
ithresh = zeros([1 stack_size(3)]);
hcmaxi = zeros([1 stack_size(3)]);
hc = zeros([stack_size(3), ops.nbins]);
for j = 1:stack_size(3)
    thisframe = img_in(:,:,j);
    hc(j,:) = histc(thisframe(:),ibins);
    [~, hcmaxi(j)] = max(hc(j,:),[],2);
    ithresh(j) = median(thisframe(thisframe>ibins(hcmaxi(j)+3)));
    if ~isempty(ops.floor)
        if ithresh(j) < ops.floor
            ithresh(j) = ops.floor;
        end
    end
end
clear thisframe ibins hc hcmaxi;

%% apply thresholds and deadframes, deadperim
for j = ops.deadframes(1)+1:stack_size(3)-ops.deadframes(2);
    mask(:,:,j) = img_in(:,:,j) > ithresh(j);
    mask(:,:,j) = imfilter(mask(:,:,j),ops.filter);
    mask(:,:,j) = mask(:,:,j)>0;
    mask(:,1:ops.deadperim(1),j) = 0;
    mask(:,stack_size(2):-1:stack_size(2)-ops.deadperim(2),j) = 0;
    mask(1:ops.deadperim(3),:,j) = 0;
    mask(stack_size(1):-1:stack_size(1)-ops.deadperim(4),:,j) = 0;
end

end