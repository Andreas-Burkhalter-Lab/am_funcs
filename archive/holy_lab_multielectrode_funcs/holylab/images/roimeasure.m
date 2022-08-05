function intens = roimeasure(im,roi1,roi2,options)
% ROIMEASURE: measure ROI intensity in an image or sequence of images
%
% Syntax:
%   intens = roimeasure(im,roi)
% where
%   im is the input image (or vimage, will be evaluated);
%   roi is a ROI structure;
% and
%   intens is a vector of mean intensities inside each ROI.
%
% If instead im is a cell array images (or vimages), then intens will be
% a nimages-by-nrois matrix, containing the intensity in each image/ROI
% combination.
%
% Alternatively, the syntax
%   intens = roimeasure(im,roi1,roi2)
% where im is a cell array of images (or vimages) will measure
% intensities in the first frame using roi1, the last frame using roi2,
% and frames in the middle with the ROI positions linearly interpolated
% between roi1 and roi2.
%
% The syntax
%   intens = roimeasure(im,roi1,roi2,'skiplast')
% will not evaluate the intensities for the last image; use this syntax
% if you're in the middle of a loop over ROI positions, to avoid
% computing twice on boundary frames.
%
% Finally, the syntax
%   intens = roimeasure(im,roi1,roi2,options)
% allows more fine-grained control through the following fields:
%   skiplast (default false): if true, will not evaluate the intensities
%     for the last image (see above);
%   progress_data: if supplied, it is a 2-vector: [offset total], where
%     total is the total number of frames being processed, and offset is
%     the frame number of the first frame supplied on this call to
%     roimeasure.
%   gaussian(default is true): if ture, use gaussian weight (centered at
%     roi, spread as sigmal =r) to caculate the integrated roi intensity
%     from an extended 8r-by-8r square region. If false, caculate 
%     the average intensity exactly from the roi region.
%
% See also: ROIPLOT, ROISTRUCT.

% Copyright 2005 by Timothy E. Holy

interpflag = 0;
if (nargin > 2 && ~isempty(roi2) )
    interpflag = 1;
end

if ~iscell(im)
    im = {im};
end
nimages = length(im);

lastimg = nimages;
if (nargin < 4)
    options = struct;
end
if (ischar(options) && strcmp(options,'skiplast'))
    options = struct('skiplast',1);
end
if ~isfield(options,'skiplast')
    options.skiplast = 0;
end
if ~isfield(options,'show_progress')
    if isfield(options,'progress_data')
        options.show_progress = 1;
    else
        options.show_progress = 0;
    end
end
if ~isfield(options,'gaussian')
  options.gaussian = 1;
end
if options.skiplast
    lastimg = nimages - 1;
end

nroi = length(roi1.type);

if (lastimg == 0)
    intens = zeros(0,nroi);
    return;
end

sz = size(im{1});

% Set up boxes containing each ROI---do this only once for the sequence
% of images. That's because sub2ind is a significant part of the
% computation.
imindx{nroi} = [];
imx{nroi} = [];
imy{nroi} = [];
for i = 1:nroi
    switch roi1.type(i)
        case 'c'
          if options.gaussian
            rfactor = 4;
          else
            rfactor = 1;
          end
            [xrange,yrange] = roim_cranges(roi1, i, rfactor);
            if interpflag
                % If interpolating, make sure box covers the range throughout
                [xrange2,yrange2] = roim_cranges(roi2, i);
                xrange(1) = min(xrange(1),xrange2(1));
                yrange(1) = min(yrange(1),yrange2(1));
                xrange(2) = max(xrange(2),xrange2(2));
                yrange(2) = max(yrange(2),yrange2(2));
            end

            if(xrange(1)<1) xrange(1)=1; end
            if(yrange(1)<1) yrange(1)=1; end
            if(xrange(2)>sz(2)) xrange(2)=sz(2); end
            if(yrange(2)>sz(1)) yrange(2)=sz(1); end

            % Encode the position of each point within box
            x = floor(xrange(1)):ceil(xrange(2));
            y = floor(yrange(1)):ceil(yrange(2));
            xlist = repmat(x,1,length(y));   % This is of the dimensions of
            % the box, holding the x-coordinate
            ylist = repmat(y,length(x),1);   % Ditto for the y-coordinate
            imx{i} = xlist(:);
            imy{i} = ylist(:);
            % Also convert to a set of indices
            imindx{i} = sub2ind(sz,imy{i},imx{i});
            % If not interpolating, also compute the weights
            if ~interpflag
                % Compute center and radius in new coordinates; get radius by
                % also transforming (x+r,y+r)
                [x,y,r] = roim_ccoords(roi1,i);
                imw{i} = roim_cweight(x-imx{i},y-imy{i},r,options.gaussian);
                imw{i} = imw{i}/sum(imw{i});
            end
        otherwise
            error(['ROI type ' roi1.type(i) ' not implemented']);
    end
end

% Loop over frames
intens(lastimg,nroi) = 0;   % Pre-allocate
for i = 1:lastimg
    if options.show_progress
        progress_bar(struct('max',options.progress_data(2),...
            'progress',options.progress_data(1)+i-1,...
            'what','Computing intensities...'));
    end
    tmpim = im{i};
    if isa(tmpim,'vimage')
        tmpim = eval(tmpim);
    end
    % Loop over ROIs
    for j = 1:nroi
        switch roi1.type(j)
            case 'c'
                if interpflag
                    % Calculate the weights
                    [x1,y1,r] = roim_ccoords(roi1,j);
                    [x2,y2,r] = roim_ccoords(roi2,j);
                    % Interpolate the position based on time
                    x = x1*(nimages-i+1)/nimages + x2*i/nimages;
                    y = y1*(nimages-i+1)/nimages + y2*i/nimages;
                    w = roim_cweight(x-imx{j},y-imy{j},r,options.gaussian);
                    w = w/sum(w);
                else
                    w = imw{j};
                end
        end % 'otherwise' will have already been caught
        intens(i,j) = sum(w.*double(tmpim(imindx{j})));
    end  % Loop over ROIs
end  % Loop over frames
if options.show_progress
    progress_bar(struct('max',options.progress_data(2),...
        'progress',options.progress_data(2),...
        'what','Computing intensities...'));
end


function [xrange,yrange] = roim_cranges(roi,i,rfactor)
    % rfactor = 1;
    % Go out as far as 4*r because we might be looking for small changes in
    % intensity, and we don't want discontinuities at places where the roi
    % definition switches.
    xrange = roi.x(i)*[1 1] + rfactor*roi.xyradius(i)*[-1 1];
    yrange = roi.y(i)*[1 1] + rfactor*roi.xyradius(i)*[-1 1];
    % Transform to actual coordinates
    [xrange,yrange] = tformfwd(roi.tform,xrange,yrange);

function [x,y,r] = roim_ccoords(roi,i)
    x = roi.x(i)+[0 roi.xyradius(i)];
    y = roi.y(i)+[0 roi.xyradius(i)];
    [x,y] = tformfwd(roi.tform,x,y);
    r = diff(x);
    x = x(1);
    y = y(1);

function w = roim_cweight(x,y,r,gaussian)
    if gaussian
      w = exp(-(x.^2 + y.^2)/(2*(0.25*r)^2))/(4*pi*r^2);
    else
      w = ~(x.^2 + y.^2 > r.^2);
    end
    
    
    
