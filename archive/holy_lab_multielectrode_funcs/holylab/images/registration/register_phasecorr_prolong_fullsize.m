function uh = register_phasecorr_prolong_fullsize(u,options)
% REGISTER_PHASECORR_PROLONG_FULLSIZE: bring u to full size of image
%
% Syntax:
%   uh = register_phasecorr_prolong_fullsize(u,options)
% where
%   u is the deformation (see register_phasecorr_improve or
%     register_phasecorr_refine);
%   options comes from register_phasecorr_initialize. You can add one extra
%     field, 'shift', to supply an additional translation to the entire
%     image.  The units of the shift should be pixels, the same as for u.
% and
%   uh is the deformation, interpolated to the full size of the image.
%
% See also: register_phasecorr_warp.

% Copyright 2011 by Timothy E. Holy

  shifting = isfield(options,'shift');

  %% Interpolate u up to full size
  uh = cell(1,options.n_dims);
  if isempty(u)
    szu = 0;
  else
    szu = size(u{1});
  end
  
  if isequal(szu,options.sz_spatial)
    uh = u;
    if shifting
      for dimIndex = 1:options.n_dims
        uh{dimIndex} = uh{dimIndex} + options.shift(dimIndex);
      end
    end
  else
    %% u is defined on a coarser grid than the image, so fill it out
    if all(szu == 1) || all(szu == 0)
      %% This is a global shift (life is easy!)
      % Check for zero shift (if so, just return original image)
      if ~all(szu == 0)
        shift = cat(2,u{:});
      else
        shift = zeros(1,options.n_dims);
      end
      if shifting
        shift = shift + options.shift;
      end
      if all(shift == 0)
        uh = [];
        return
      end
      % Nonzero shift, so replicate the shift to the full size of the image
      for dimIndex = 1:options.n_dims
        uh{dimIndex} = repmat(shift(dimIndex),options.sz_spatial);
      end
    else
      %% This is a deformation, interpolate
      if isfield(options,'pyramid')
        % We can use array_prolong, which is much faster
        matched = false;
        for sizeIndex = 2:length(options.pyramid)
          if all(szu == options.pyramid(sizeIndex).sz)
            matched = true;
            break
          end
        end
        if ~matched
          error('Using a pyramid, but failed to match the size.');
        end
        for dimIndex = 1:options.n_dims
          utmp = u{dimIndex};
          if shifting
            utmp = utmp + options.shift(dimIndex);
          end
          for levelIndex = sizeIndex-1:-1:1
            utmp = array_prolong(utmp,options.pyramid(levelIndex).sz);
          end
          uh{dimIndex} = utmp;
        end
      else
        % Generate by interpolation
        ux = cell(1,options.n_dims);
        for dimIndex = 1:options.n_dims
          n = szu(dimIndex);
          N = options.sz_spatial(dimIndex);
          % % Define u on an evenly-spaced grid inset within the image
          % ux{dimIndex} = cast((0.5:N-0.5)/N * (n+1),options.class);
          % Define u on an evenly-spaced grid that includes the corners
          ux{dimIndex} = cast((0:N-1)*((n-1)/(N-1))+1,options.class);
        end
        [ux{:}] = ndgrid(ux{:});
        % Calculate the full-sized u grid
        for dimIndex = 1:options.n_dims
          utmp = u{dimIndex};
          if shifting
            utmp = utmp + options.shift(dimIndex);
          end
          uh{dimIndex} = iminterp(utmp,ux{:},'extrap');
        end
      end
    end
  end
  