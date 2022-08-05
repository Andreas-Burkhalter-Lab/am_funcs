classdef opticalRayBundle
% opticalRayBundle: a collection of rays in a common coordinate system
%
% This defines rays via a position, direction of propagation, wavelength,
% and intensity.
%
% Constructor:
%   rb = opticalRayBundle(p,e,wavelength)
%   rb = opticalRayBundle(p,e,wavelength,intensity)
% Defines the rays by their starting position (p, a 3-by-N matrix), their
% direction of propagation (e, a 3-by-N matrix, does not have to be
% normalized), and wavelength in vacuum (either a scalar or a 1-by-N
% vector, in units of nanometers).
% Optionally, you can specify an intensity for each ray (a 1-by-N
% vector).
% See also optical_collimated_bundle and optical_point_bundle for
% alternative ways to construct a ray bundle.
%
% Methods:
%   rb.valid = validFlag;
% Set the "valid flag" for the rays. This is handy for keeping track of
% rays that fail to make it through an aperture, for example. It subsets
% from the currently-valid set, so you only have to concern yourself with
% "turning off" rays as they become invalid.
%
%   p = rb.p;
% Obtain the current location, just for valid rays.
%   p = rb.p_all;
% Obtain the current location of all rays (even invalid ones). This is
% mostly useful for contatenation purposes, when you need the size of the
% outputs to be fixed.
%   e = rb.e;
%   e = rb.e_all;
% Obtain the current propagation direction (first one returns data just for
% valid rays). e is guaranteed to be normalized.
%   L = rb.opl;
% Obtain the total optical path length traversed by each valid ray.
% And so forth (see available properties).
%
%   [rb,p] = rb.propagate(t,n);
% Propagate each ray a distance t (a 1-by-N_valid vector) along its
% propagation direction, through a medium of refractive index n. n can be a
% scalar (applied to all rays), or a vector of length N_valid, or a vector
% of length N_total.
% The optional output p contains the new positions of all valid rays.
%
%   rb = rb.refract(normal,nratio,tirFlag);
% Refract each ray through a surface specified by the normal and refractive
% index ratio nratio = npre./npost between the input and output sides,
% respectively. normal is a 3-by-N_valid matrix (will be normalized for
% you, so you don't have to worry about normalization). If tirFlag is true,
% then rays that experience total internal reflection will "become" the
% reflected ray; if tirFlag is false (the default), such rays are set to be
% invalid.
% Note that "normal" needs to be supplied in the same coordinate system as
% the rays ("absolute" coordinates), which may not be the same as the
% "relative" coordinates that conveniently assume the optic axis is along
% [0 0 1]. Any component that inherits from opticalAperature should use the
% "transform_rotate_inverse" method to convert normal to absolute
% coordinates.
% 
%   [mean_coefs,cov_coefs] = rb.mean_covariance(normal)
% Calculate the intensity-weighted mean and covariance of ray position 
% at intersection with planes perpendicular to 'normal'. The mean is
% parametrized as
%   m = z*mean_coefs(:,1) + mean_coefs(:,2),
% where mean_coefs is 3-by-2, and the covariance as
%   C = z^2*cov_coefs(:,:,1) + z*cov_coefs(:,:,2) + cov_coefs(:,:,3).
% In these expressions, z is the normal displacement from the base plane,
% defined as the plane that contains the z=0 point mean_coefs(:,2).
%   For beams that are not collimated, the information provided by this
% function can be used to find the "center of focus." The value of z that
% minimizes total radius is
%      zmin = - trace(cov_coefs(:,:,2))/(2*trace(cov_coefs(:,:,1))
% Alternatively, if you want to minimize the width along just the axis
% defined by the column vector v, then
%      zmin = - (v'*cov_coefs(:,:,2)*v)/(2*v'*cov_coefs(:,:,1)*v).
% Naturally, C will have a zero eigenvalue along the axis 'normal'.
%
% See also: optical_collimated_bundle, optical_point_bundle.

% Copyright 2010 by Timothy E. Holy

  properties (SetAccess = private,GetAccess = public)
    % Properties of the rays themselves
    p_all = [];
    e_all = [];
    wavelength_all = [];
    intensity_all = [];
    opl_all = [];
  end
  properties (Dependent = true, GetAccess = public, SetAccess = private)
    % These just return values for the valid rays
    p
    e
    wavelength
    intensity
    opl
  end
  properties (Access = private)    
    validFlag = [];
  end
  properties (Dependent = true, GetAccess = private, SetAccess = public)
    valid  % for setting validFlag; subsets from the previously-valid rays
  end
  methods
    %% Constructor
    function rb = opticalRayBundle(pin,ein,w,Iin)
      if (nargin > 0)
        [d,N] = size(pin);
        if (d ~= 3)
          error('Position matrix p must be 3-by-N');
        end
        rb.p_all = pin;
        if ~isequal(size(ein),[d N])
          error('e must have the same size as p');
        end
        enorm = sum(ein.^2,1);
        if any(enorm == 0)
          error('Propagation direction cannot have zero length');
        end
        rb.e_all = bsxfun(@rdivide,ein,sqrt(enorm));
        if isscalar(w)
          rb.wavelength_all = repmat(w,1,N);
        elseif (numel(w) == N)
          rb.wavelength_all = w(:)';
        else
          error('wavelength must be a scalar or 1-by-N');
        end
        if (nargin > 3)
          if isscalar(Iin)
            rb.intensity_all = repmat(Iin,1,N);
          elseif (numel(Iin) == N)
            rb.intensity_all = Iin(:)';
          else
            error('Intensity must be a scalar or 1-by-N');
          end
        else
          rb.intensity_all = ones(1,N);
        end
        rb.opl_all = zeros(1,N);
        rb.validFlag = true(1,N);
      end
    end  % end constructor
    
    %% p
    function p1 = get.p(rb)
      p1 = rb.p_all(:,rb.validFlag);
    end
    
    %% e
    function e1 = get.e(rb)
      e1 = rb.e_all(:,rb.validFlag);
    end
    
    %% wavelength
    function w = get.wavelength(rb)
      w = rb.wavelength_all(rb.validFlag);
    end
    
    %% intensity
    function i = get.intensity(rb)
      i = rb.intensity_all(rb.validFlag);
    end
    
    %% optical path length
    function Lo = get.opl(rb)
      Lo = rb.opl_all(rb.validFlag);
    end
    
    %% valid
    function rb = set.valid(rb,validF)
      validIndex = find(rb.validFlag);
      rb.validFlag(validIndex(~validF)) = false;
    end
    
    %% propagate
    function [rb,pout] = propagate(rb,t,n)
      n_valid = sum(rb.validFlag);
      if (numel(t) ~= n_valid)
        error('The size of t does not agree with the number of valid rays');
      end
      if isscalar(n)
        n = repmat(n,1,n_valid);
      elseif (numel(n) == length(rb.validFlag))
        % n supplied for all rays, but we just want valid ones
        n = n(rb.validFlag);
      end
      if (numel(n) ~= n_valid)
        error('The size of n does not agree with the number of valid rays');
      end
      rb.opl_all(:,rb.validFlag) = rb.opl + t.*n;
      pout = rb.p + bsxfun(@times,t,rb.e);
      rb.p_all(:,rb.validFlag) = pout;
    end  % end propagate
    
    %% refract
    function rb = refract(rb,normal,nratio,tirFlag)
      if (nargin < 4)
        tirFlag = false;
      end
      n_valid = sum(rb.validFlag);
      if (size(normal,2) ~= n_valid)
        error(['The number of normal vectors must be equal to the number' ...
          ' of valid rays']);
      end
      if (numel(nratio) == length(rb.validFlag))
        % nratio was supplied for all rays, but we just want valid ones
        nratio = nratio(rb.validFlag);
      end
      if ~(isscalar(nratio) || numel(nratio) == n_valid)
        error(['The number of index ratios must be a scalar or equal to' ...
          ' the number of rays']);
      end
      % Normalize the normal
      normnorm = sum(normal.^2,1);
      if any(normnorm == 0)
        error('None of the normals may have zero length');
      end
      normal = bsxfun(@rdivide,normal,sqrt(normnorm));
      % Compute the dot product with the normal
      eout = rb.e;
      c1 = -sum(normal.*eout,1);
      s = sign(c1);
      % Check for total internal reflection
      c2arg = 1 - nratio.^2.*(1-c1.^2);
      refractFlag = c2arg >= 0;   % if false, TIR
      if ~all(refractFlag)
        if tirFlag
          % For rays with TIR, replace them with their reflected version
          nrF = ~refractFlag;
          eout(:,nrF) = eout(:,nrF) + 2*bsxfun(@times,c1(nrF),normal(:,nrF));
        else
          % Invalidate the TIR rays
          rb.valid = refractFlag;
        end
      end
      % Snell's law
      if ~isscalar(nratio)
        nratio = nratio(refractFlag);
      end
      c2 = sqrt(c2arg(refractFlag));
      eout = bsxfun(@times,nratio,eout(:,refractFlag)) + ...
        bsxfun(@times,nratio.*c1(refractFlag) - c2.*s(refractFlag),normal(:,refractFlag));
      eout = bsxfun(@rdivide,eout,sqrt(sum(eout.^2,1)));
      rb.e_all(:,rb.validFlag) = eout;
    end % end refract

    %% mean_covariance
    function [mean_coefs,cov_coefs] = mean_covariance(rb,normal)
      p = rb.p;
      e = rb.e;
      I = rb.intensity;
      normal = normal(:);
      edotn = sum(bsxfun(@times,e,normal),1);
      if any(edotn == 0)
        error('At least one ray is perpendicular to "normal," so it is not possible to provide statistics along this axis');
      end
      zcoef = bsxfun(@rdivide,e,edotn);  % coefficient of ray position as fcn of z
      Isum = sum(I);
      % Calculate the mean of the ray base points
      pmean = sum(bsxfun(@times,I,p),2)/Isum;
      % Calculate the component of displacement from pmean in the plane
      % containing pmean
      dp = bsxfun(@minus,p,pmean);
      t = -sum(bsxfun(@times,dp,normal),1)./edotn;
      dp = dp + bsxfun(@times,t,e); % "propagate" displacement to that plane
      % Calculate the slope of mean position with z
      zcoefmean = sum(bsxfun(@times,I,zcoef),2)/Isum;
      mean_coefs = [zcoefmean pmean];
      % Calculate the covariance
      dzcoef = bsxfun(@minus,zcoef,zcoefmean);
      cov_coefs = cat(3,bsxfun(@times,I,dzcoef)*dzcoef',...
        bsxfun(@times,2*I,dp)*dzcoef',...
        bsxfun(@times,I,dp)*dp')/Isum;
    end    
  end % methods
end
