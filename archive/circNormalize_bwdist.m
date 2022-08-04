%% this function is an attempt at speeding up with bwdist... but is slower

%%% divide each pixel value by the mean value of pixels within a
%
% normedImage = circNormalize(im,umRadius,zoom,scope,savetifPars)
%
%%% pixRadius-size circle around it
% updated 6/13/18

function normedImage = circNormalize_bwdist(im,umRadius,zoom,scope,savetifPars)


if ~exist('scope','var') || isempty(scope) %  1392*1040
    scope = 'epimicro';
end

pixRadius = round(umRadius .* pixPerUm(zoom,scope));
 
 
if ischar(im)
    imfile = im;
    im = imread(im);
else
    imfile = '';
end
im = im(:,:,1);

% draw circle for determining which pix are within pixRadius
boxlength = 2*sqrt(2)*pixRadius;
[xmesh ymesh] = meshgrid(1:boxlength,1:boxlength);
meshcntr = (boxlength+1)/2;
vals = sqrt((xmesh-meshcntr).^2 + (ymesh-meshcntr).^2); % for drawing circle
circImage = vals <= pixRadius; % 1s are in the circle, 0s outside
[circsubs(:,1) circsubs(:,2)] = ind2sub(size(circImage),find(circImage)); % circle as subscripts
circsubs = circsubs - round(meshcntr); % center at zero

im = double(im);
dims = size(im);
normedImage = NaN(size(im));
npix = numel(im);
% % % origmax = max(max(im));
origmean = mean(mean(im));

barhandle = waitbar(0,['Normalizing ' imfile '...']);
onepointim = false(size(im));
for indpix = 1:npix
    if indpix>1
        onepointim(indpix-1) = false;
    end
    onepointim(indpix) = true;
    circmean = mean(im(bwdist(onepointim) > pixRadius)); 
    normedImage(indpix) = im(indpix) / circmean;

    if rem(indpix,1e5) == 0
        try waitbar(indpix/npix,barhandle); end
    end
end
try close(barhandle); end

% scale resulting values to be have same mean value as original
% % % % rescale_factor = origmax / max(max(normedImage)); %%%% changed to mean 12/7/17  
naninds = find(isnan(normedImage));
normedImage(naninds)  = 0;
rescale_factor = origmean / mean(mean(normedImage));
normedImage = normedImage .* rescale_factor;

if ~exist('savetifPars','var')
    savetifPars = struct;
end

if ~isempty(imfile)
    if ~isempty(fileparts(imfile)) % if includes dir
        savetif(normedImage, [fileparts(imfile) filesep getfname(imfile) '_normed' num2str(pixRadius) 'pix' num2str(umRadius) 'um'],savetifPars)        
    else
        savetif(normedImage, [getfname(imfile) '_normed' num2str(pixRadius) 'pix' num2str(umRadius) 'um'],savetifPars)
    end
end
