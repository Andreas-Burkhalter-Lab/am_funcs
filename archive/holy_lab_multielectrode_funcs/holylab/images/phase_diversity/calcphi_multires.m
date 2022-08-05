function phi = calcphi_multires(im,H0,phi0,options)
% CALCPHI_MULTIRES: use multiresolution strategy for phase diversity
%
% Takes advantage of noise-reduction and faster computation by starting
% from reduced-size images and then work to fine-scale images.
%
% Syntax:
%   phi = calcphi_multires(im,H0,phi0,options)
% where the syntax is identical to CALCPHI2D, with the restriction that
% the parameters are actual phi matrices (e.g., use CALCPHI_LINEARSUMS
% rather than CALCPHI_ZERNIKE.)  "options" may have the following extra
% fields:
%   quit_early (default false): don't actually use the finest-scale
%     resolution, merely interpolate up to this scale
% 
% One enhancement of this algorithm is that your images do not have to be
% the same size as H0 and phi0---it will take subregions of the images to
% adjust appropriately. (This permits one to handle edge effects by
% focusing on the middle of the data. Edge effects are still important,
% but may be (?) less critical.)
%
% See also: CALCPHI2D.

% Copyright 2009 by Timothy E. Holy

  if (nargin < 4)
    options = struct;
  end
  options = default(options,'quit_early',false);
  
  %% Multiresolution data preparation
  data = cell(0,3);
  Hs = fftshift(H0);
  phi0s = [];
  if ~isempty(phi0)
    phi0s = fftshift(phi0);
  end
  isdone = false;
  while ~isdone
    % Process everything related to the pupil
    % Under a change in resolution, the pupil simply "shaves off" the
    % blocked region
    data{end+1,2} = fftshift(Hs);
    data{end,3} = fftshift(phi0s);
    sz = size(Hs);
    newsz = ceil((sz+1)/2);
    szdiff = sz - newsz;
    gap = round(szdiff/2);
    indx = {gap(1)+(1:newsz(1)), gap(2)+(1:newsz(2)), ':'};
    Hs = Hs(indx{:});
    if ~isempty(phi0s)
      phi0s = phi0s(indx{:});
    end
    % Process the image. If the images are larger than the pupil, cut out
    % the middle region (but save everything so that we handle edges well,
    % or at least try to!).
    imsz = size(im);
    imgap = round((imsz(1:2) - sz)/2);
    indx = {imgap(1)+(1:sz(1)), imgap(2)+(1:sz(2))};
    data{end,1} = im(indx{:},:);
    im = array_restrict(im,[true true false]); % smaller for next time
    % Check to see if the edge of the pupil is intersecting the
    % boundaries of the image. If so, we shouldn't continue.
    isdone = any(Hs(:,1) ~= 0) || any(Hs(1,:) ~= 0);
  end
  

  %% Multiresolution phi optimization
  n_levels = size(data,1);
  for level = n_levels:-1:1
    if (options.quit_early && level == 1 && n_levels > 1)
      % Quit early!
      phi = data{1,3};
      return
    end
    phi = calcphi2d(data{level,:},options);  % Do the optimization at the current level
    if (level > 1)
      % Insert the result as the starting guess for the next level of
      % resolution
      sz = size(phi);
      extra_sz = [];
      if (length(sz) > 2)
        extra_sz = sz(3:end);
        sz = sz(1:2);
      end
      sznext = size(data{level-1,2});  % if we use H don't have to worry empty
      phis = fftshift(phi);
      gap = round((sznext - sz)/2);
      phisnew = zeros([sznext extra_sz]);
      indx = {gap(1)+(1:sz(1)),gap(2)+(1:sz(2)),':'};
      phisnew(indx{:}) = phis;
      data{level-1,3} = fftshift(phisnew);
    end
  end
end
