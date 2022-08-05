function params_out = tune_dm2d(params_in,rays,target_angles)
% params: a structure array with the fields described in
%   OPT2DDMFLAT. Currently you can supply 1 (goal: 2) elements in this
%   structure array, allowing up to 2 (sequential) DMs.
% rays: a structure array describing rays that are just about to strike
%   the first DM. See RAY.
% target_angles: the desired angles (with the horizontal) that each ray
%   would ideally have. These angles are specified in terms of the rays
%   _before_ striking the DM, i.e., what they would be in the absence of
%   the DM(s). The DMs might introduce an overall reflection angle that
%   is not to be included in this estimate.
  
  n_dms = length(params_in);
  n_rays = length(rays);
  if (n_dms > 1)
    error('More than 1 DM is not currently supported');
  end

  % Trace the rays forwards to their strike position with the first DM
  th = params_in(1).theta;
  ptmp = struct('center',params_in(1).center,'normal',[-sin(th);cos(th)]);
  [t,strike_pos1] = ray_flat_intersection(ptmp,rays);
  [strikeIndex1,strikeFrac1] = tdm_strikel(params_in(1),strike_pos1);

  % Calculate each ray's input angle
  e = [rays.e];
  cosra = e(1,:);
  nzFlag = cosra ~= 0;
  ray_angle(nzFlag) = atan(e(2,nzFlag) ./ cosra);
  if ~all(nzFlag)
    % These would have given a divide by zero warning
    ray_angle(~nzFlag) = ((e(2,~nzFlag) > 0)-0.5) * pi;
  end
  
  if (n_dms == 1)
    % Simple case: just calculate the phis directly as
    %     phi(rayi) = -delta(rayi)/2
    % The one subtlety is that phi(rayi) is determined by a finite set of
    % parameters, and so this equation may not be solvable exactly---use
    % SVD to solve in a least-squares sense.  At any position, phi(rayi)
    % is determined by linear interpolation between the two nodal points
    % on either side.
    delta = target_angles - ray_angle;
    A = sparse([1:n_rays 1:n_rays],...
	       [strikeIndex1 strikeIndex1+1],...
	       [1-strikeFrac1 strikeFrac1],...
	       n_rays,length(params_in.DMparams));
    params_out.DMparams = -(A\delta(:))/2;
  end
  
%   theta: the angle of the "plane" of the DM with respect to horizontal
%   DMparams: the local tilts of the DM
%   center

function [strikeIndex,strikeFrac] = tdm_strikel(params_in,strike_pos)
  n_rays = size(strike_pos,2);
  parallel = [cos(params_in.theta);sin(params_in.theta)];
  P = diag(parallel);
  dx = strike_pos - repmat(params_in.center(:),1,n_rays);
  l = sum(P*dx);
  % Convert this to an index & frac
  lbreak = linspace(-params_in.length/2,params_in.length/2, ...
		    length(params_in.DMparams));
  dl = diff(lbreak(1:2));
  [n_per_segment,strikeIndex] = histc(l,lbreak);
  struckFlag = strikeIndex > 0;
  strikeFrac = nan(size(strikeIndex));
  strikeFrac(struckFlag) = ...
      (l(struckFlag) - lbreak(strikeIndex(struckFlag)))/dl;
  