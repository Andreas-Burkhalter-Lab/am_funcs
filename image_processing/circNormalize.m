%%% divide each pixel value by the mean value of pixels within a
%
% normedImage = circNormalize(im,umRadius,zoom,scope,savetifPars)
%
%%% pixRadius-size circle around it
% updated 18/06/28 on thermaltake

function normedImage = circNormalize(im,umRadius,zoom,scope,savetifPars)


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
origmean = mean(mean(im));

barhandle = waitbar(0,['Normalizing ' imfile '...']);
for indpix = 1:npix
    [thiscoord(1) thiscoord(2)] = ind2sub(dims,indpix);
    inRadiusSubs = [circsubs(:,1)+thiscoord(1), circsubs(:,2)+thiscoord(2)];
    if any(thiscoord < pixRadius+1) || any(dims-thiscoord < pixRadius+1) % if inRadiusSubs are outside of image
        deleteRows = inRadiusSubs(:,1)<1 | inRadiusSubs(:,2)<1  |  inRadiusSubs(:,1)>dims(1) | inRadiusSubs(:,2)>dims(2);
        inRadiusSubs(deleteRows,:) = [];
    end
    
    %% norming loop
    circinds = sub2ind(size(im), inRadiusSubs(:,1), inRadiusSubs(:,2));
    circvals = im(circinds);
    circmean = mean(circvals);  
    normedImage(indpix) = im(indpix) / circmean;
%%

    if rem(indpix,1e5) == 0
        try waitbar(indpix/npix,barhandle);end
    end
end
try close(barhandle); end

% scale resulting values to be have same mean value as original
naninds = find(isnan(normedImage));
normedImage(naninds)  = 0;
rescale_factor = origmean / mean(mean(normedImage));
normedImage = normedImage .* rescale_factor;

if ~exist('savetifPars','var')
    savetifPars = struct;
end
if strcmp(scope,'epimicro_old') && ~isfield(savetifPars,'conv_factor') % if image is dark and conv_factor not specified
    savetifPars.conv_factor = 16;% brighten epimicro old images, which were saved at low absolute pixel values
end

if ~isempty(imfile)
    if ~isempty(fileparts(imfile)) % if includes dir
        savetif(normedImage, [fileparts(imfile) filesep getfname(imfile) '_normed' num2str(pixRadius) 'pix' num2str(umRadius) 'um'],savetifPars)        
    else
        savetif(normedImage, [getfname(imfile) '_normed' num2str(pixRadius) 'pix' num2str(umRadius) 'um'],savetifPars)
    end
end


% % % % % %%% use filter; results in negative values
% % % % % flt = -(1/[length(find(circImage))-1]) * ones(size(circImage));
% % % % % flt(~circImage) = 0;
% % % % % flt(round(meshcntr),round(meshcntr)) = 1;
% % % % % normedImage = imfilter(im,flt,'symmetric');



% % % % %%%%%%%%%%%% slower version of norming loop
% % % %     if indpix == 1
% % % %         circinds = sub2ind(size(im), inRadiusSubs(:,1), inRadiusSubs(:,2));
% % % %         circvals = im(circinds);
% % % %         circmean = mean(circvals);
% % % %     elseif indpix >= 1
% % % %         oldinds = circinds;
% % % %         oldcircmean = circmean;
% % % %         circinds = sub2ind(size(im), inRadiusSubs(:,1), inRadiusSubs(:,2));
% % % %         if length(oldinds) ~= length(circinds) % use slow method if normalization pool sizes don't match
% % % %             circvals = im(circinds);
% % % %             circmean = mean(circvals);
% % % %         elseif length(oldinds) == length(circinds) % add new values to be used in normalizing to values shared from the last pixel
% % % %             ncircinds = length(circinds);
% % % %             [~, oldshared, newshared] = intersect(oldinds,circinds);
% % % %             oldunique = oldinds;
% % % %             oldunique(oldshared) = [];
% % % %             olduniquevals = im(oldunique);
% % % %             newunique = circinds;
% % % %             newunique(newshared) = [];
% % % %             newuniquevals = im(newunique);
% % % %             %%%% subtract old nonshared values from mean, add new unshared
% % % %             circmean = oldcircmean - sum(olduniquevals)/ncircinds + sum(newuniquevals)/ncircinds;
% % % %         end
% % % %     normedImage(indpix) = im(indpix) / circmean;
% % % %     end
