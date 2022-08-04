%%%% save matlab data image as tif
% savetif(img,filename,pars,tagstruct)
%   pars = pars struct; alternatively, a single value to use as conversion factor
%       pars options = conv_factor, keepbitdepth, zero_transparent
%%%%%%%%%%% note: to use alpha channel transparency in illustrator, save image as .psd then open in illustrator
%%%%%%% updated 2020/10/11 on thermaltake

function savetif(img,filename,pars,tagstruct)

imgclass = class(img);
tagstruct = vardefault('tagstruct',struct);
pars = vardefault('pars',struct);
if isnumeric(pars)
    c = pars;
    pars = struct;
    pars.conv_factor = c;
elseif ~isstruct(pars)
    error('bad third arg')
end
    
if islogical(img) || max(img(:)) < 10 %%% if pixel values are logical or very low
    pars.conv_factor = field_default(pars,'conv_factor', 16^3.3); % brighten white spots to make visible in photoshop; 16^3.3 good for maxval=6
    pars.zero_transparent = field_default(pars,'zero_transparent',true);
else
%     pars.conv_factor = field_default(pars,'conv_factor', 16); % for epimicro old images, which were saved at low absolute pixel values
    pars.conv_factor = field_default(pars,'conv_factor', 1); 
end


if ~isfield(tagstruct,'SamplesPerPixel')
    tagstruct.SamplesPerPixel = size(img,3);
end

% set image to rgb or grayscale
if ~isfield(tagstruct,'Photometric') && size(img,3)==3 && tagstruct.SamplesPerPixel==3
    tagstruct.Photometric = Tiff.Photometric.RGB;
else
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
end

% decide whether an alpha layer will be added to make zeros transparent
%%%% note: nonzero values will have equal alpha but it will likely NOT be max alpha, so they'll be partly transparent
if isfield(pars,'zero_transparent') && pars.zero_transparent
    isNotZero = logical(max(img,[],3)); % set as nonzero all pixels which are nonzero at any value in 3rd dimension (RGB if 3D)
    alphamask = max(img,[],3);
    alphamask(isNotZero) = max(img(:)); % set nonzeros to equal alpha
    img(:,:,size(img,3)+1) = alphamask; % add alpha layer
    tagstruct.SamplesPerPixel = tagstruct.SamplesPerPixel + 1; % indicate there's an extra alpha layer
    tagstruct.ExtraSamples = Tiff.ExtraSamples.AssociatedAlpha; 
end

tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.ImageLength = size(img,1);
tagstruct.ImageWidth = size(img,2);

img = pars.conv_factor * double(img);

if ~exist('pars','var') || ~isfield(pars,'keepbitdepth' ) || ~pars.keepbitdepth %%% convert to 16 bit depth
    tagstruct.BitsPerSample = 16;
    img = uint16(img);
else %%%% keep original bit depth
    switch imgclass
        case 'uint8'
            tagstruct.BitsPerSample = 8;
            img = uint8(img);
        case 'uint16'
            tagstruct.BitsPerSample = 16;
            img = uint16(img);
        otherwise
            error('unrecognized image bit depth')
    end
end

if ~exist('filename','var')
    filename = ['newtif_' date '_' num2str(now)];
end

if length(filename)>3 && strcmp(filename(end-3:end),'.tif')
    filename(end-3:end) = [];
end
t = Tiff([filename '.tif'],'w');
t.setTag(tagstruct)
% % % t.setTag('ExtraSamples',Tiff.ExtraSamples.AssociatedAlpha);
t.write(img);