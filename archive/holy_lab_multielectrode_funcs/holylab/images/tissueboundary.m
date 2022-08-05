function [z,zi] = tissueboundary(stk,params)
% TISSUEBOUNDARY: find the top surface of the tissue in fluorescence
% image stacks
%
% This function calculates the surface from a set of grid points,
% interpolating the height at points between the grid points.  There is
% then an energy functional
%     E = (pixel intensity penalty) + (potential energy)
% where the pixel intensity penalty sums up the total pixel intensity
% above the current surface, and the potential energy term is
% proporational to the integrated height of the surface.  Thus, energy is
% decreased by lowering the surface, at least until one starts paying a
% significant price in terms of pixel intensity.
%
% Syntax:
%   pixe = tissueboundary(stk)
%   z = tissueboundary(stk,params)
%   [z,zi] = tissueboundary(...)
% where
%   stk is an image stack, with the tissue oriented so that the surface
%     is approximately perpendicular to the first coordinate, and
%     moving above the tissue increases the value of the first
%     coordinate;
%   params is a structure which may have the following fields:
%     height2intens: the tradeoff between pixel intensity and potential
%       energy.  You'll have to experiment to figure out how to set this
%       value.  Larger values mean that the algorithm will be less
%       inclined to permit bright patches above the surface.
%       If you omit this field, then the algorithm will return the summed
%       pixel intensity vs. height.  You can then use this information to
%       help you set height2intens so that the minimum energy occurs at the
%       boundary.
%     z0: if supplied, used as the starting guess & grid for z (so gridsep
%       is ignored)
%     zmin (default 1): a value safely below the tissue
%       surface---you can use this to insure that the surface does not
%       end up below the tissue during the global-scalar phase
%     gridsep(default 50): the spacing (in pixels) between grid points,
%       the same in x & y
%     interactive (default true): if true, plots the energy and waits for
%       user confirmation of parameters
%     plot (default false): if true, plots the surface at each iteration
%     backoff (default 0): the amount to increase z by after finding the
%       optimal scalar value, before going to the grid parametrization of
%       z
%     maxiter: if present, sets the maximum number of iterations of the
%       minimization algorithm.  Otherwise, the default value chosen by
%       Matlab is taken.
% and
%   pixe(z) is the pixel intensity at height z, averaged across the x,y
%     plane;
%   z is the value of z on the grid;
%   zi is the interpolated z for each x,y.
%
% See also: TISSUEBOUNDARY_CHECK.

% Copyright 2006 by Timothy E. Holy

% Change history:
%    5/5/2006: Big improvement in efficiency by pre-computing the
%      pixel-intensity energy---then can just lookup

% Add starting guess?  
  
  if (nargin < 2)
    params = struct;
  end
  if ~isfield(params,'zmin')
    params.zmin = 1;
  end
  if ~isfield(params,'gridsep')
    params.gridsep = 50;
  end
  if ~isfield(params,'plot')
    params.plot = false;
  end
  if ~isfield(params,'backoff')
    params.backoff = 0;
  end
  if ~isfield(params,'interactive')
    params.interactive = true;
  end
  
  sz = size(stk);
  if ~isfield(params,'height2intens')
    % Return average pixel intensity vs. height, to help user set
    % height2intens
    stk_avg = stk;
    for i = 2:ndims(stk)
      stk_avg = mean(stk_avg,i);  % This works for any dimensionality stack
    end
    z = squeeze(stk_avg);
    return
  end
  % Find the optimal value for a scalar z
  % We can compute the energy vs height efficiently using cumsum
  colons = {':'};
  colons = repmat(colons,[1 length(sz)-1]);
  stkcum = flipdim(cumsum(single(flipdim(stk(params.zmin:end,colons{:}),1)),1),1);
  pixenergy = stkcum;
  for i = 2:ndims(stk)
    %pixenergy = mean(pixenergy,i);
    pixenergy = max(pixenergy,[],i);
  end
  energy = pixenergy + params.height2intens*(1:size(stk,1))';
  energy(1:params.zmin) = Inf; % make sure these are not chosen
  if params.interactive
    hfig = figure;
    plot(energy)
    ylabel('Energy')
    xlabel('z (1st coordinate)')
    fprintf(['Make sure there is a minimum in the energy near the ' ...
             'boundary.\nPress any key to continue.\n']);
    pause
    delete(hfig)
  end
  if ~isfield(params,'z0')
    [minenergy,z] = min(energy);
    % If applicable, back z off a bit
    z = z + params.backoff;
    % Insert more landmark points and re-optimize
    zrep = ceil(sz(2:end)./params.gridsep);
    if (length(zrep) == 1)
      % Treat the case of a 2-d image specially here
      zrep = [zrep 1];
    end
    z = repmat(z,zrep);
  else
    z = params.z0 - (params.zmin-1);
  end
  M = tb_interp(size(z),sz);
  tb_energy(z,stkcum,M,params)
  % Compute optimum
  minops = optimset('display','iter','LargeScale','off','TolX',0.001,...
    'DiffMinChange',0.1,'DiffMaxChange',1);
  if isfield(params,'maxiter')
    minops.MaxIter = params.maxiter;
  end
  lb = zeros(size(z)) + params.zmin + 0.1;
  ub = zeros(size(z)) + sz(1);
  %z = fmincon(@(zc) tb_energy(zc,stkcum,M,params),z,...
  %            [],[],[],[],lb,ub,[],minops);
  z = fminsearch(@(zc) tb_energy(zc,stkcum,M,params),z,minops);
  %for i = 1:5;
  %  z = fmin1d(z,stkcum,M,params);
  %end
  zi = M*z(:);
  if ~isvector(z)
    zi = reshape(zi,sz(2:end));
  else
    zi = reshape(zi,[1 length(zi)]);
  end
  z = z + params.zmin-1;
  zi = zi + params.zmin-1;
  
function M = tb_interp(sz_z,sz_stk,params)
  % Calculate the (x,y) coordinates of the specified points (in pixel units)
  x = cell(1,length(sz_stk)-1);
  xin = x;
  for i = 1:length(x)
    xin{i} = round(linspace(1,sz_stk(i+1),sz_z(i)));
    x{i} = 1:sz_stk(i+1);
  end
  M = interpmatrix(xin{:},x{:});

function energy = tb_energy(z,stkcum,M,params)
  sz_stk = size(stkcum);
  if any(z(:) > sz_stk(1) | z(:) < params.zmin)
    energy = Inf;
    return;
  end
  colons = {':'};
  colons = repmat(colons,[1 ndims(stkcum)-1]);
  % Calculate the energy.
  if (numel(z)== 1)
    % Special case for when height is a scalar, to make things faster
    cz = ceil(z);
    fz = floor(z);
    if (cz == fz)
      energy = stkcum(z,colons{:});
    else
      % Linear interpolation for fractional z
      energy = (cz-z)*stkcum(fz,colons{:}) + (z-fz)*stkcum(cz,colons{:});
    end
    energy = sum(energy(:)) + params.height2intens*numel(energy)*z;
  else
    % From the supplied landmark points, calculate the tissue height at all
    % locations
    zi = M*z(:);
    if ~isvector(z)
      zi = reshape(zi,sz_stk(2:end));
    else
      zi = reshape(zi,[1 length(zi)]);
    end
    % Plot, if desired
    if params.plot
      if isvector(z)
        plot(z)
      else
        hs = surf(z);
        %set(hs,'EdgeColor','none')
      end
      drawnow
    end
    zif = floor(zi);
    % Use indexing efficiently: calculate the index corresponding to each
    % (z,x,y) point on the boundary
    x = cell(1,ndims(stkcum)-1);
    X = x;
    for i = 1:ndims(stkcum)-1
      x{i} = 1:size(stkcum,i+1);
    end
    if (length(x) > 1)
      [X{:}] = ndgrid(x{:});
    else
      X = x;
    end
    boundaryindex = sub2ind(size(stkcum),zif,X{:});
    energy = (zif(:)+1 - zi(:)).*stkcum(boundaryindex(:)) + ...
             (zi(:) - zif(:)).*stkcum(min(boundaryindex(:)+1,numel(stkcum)));
    energy = sum(energy) + params.height2intens*sum(zi(:));
  end
  
  
function z = fmin1d(z0,stkcum,M,params)
  z = z0;
  for i = 1:numel(z0)
    z(i) = fminbnd(@(x) fmin1d_inner(x,z0,i,stkcum,M,params),1,size(stkcum,1));
  end
  
function energy = fmin1d_inner(x,z,i,stkcum,M,params)
  z(i) = x;
  energy = tb_energy(z,stkcum,M,params);