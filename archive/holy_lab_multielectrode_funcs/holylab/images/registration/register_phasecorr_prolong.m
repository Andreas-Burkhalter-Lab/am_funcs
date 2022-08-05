function [sz_unew,uh] = register_phasecorr_prolong(u,options)
% register_phasecorr_prolong: resize to a finer grid
%
% Syntax:
%   sz_unew = register_phasecorr_prolong(u,options)
%   [sz_unew,uh] = register_phasecorr_prolong(u,options)
% where
%   u is a deformation, expressed either as a cell array (one entry for each
%     spatial coordinate of the image) or an array (with the spatial
%     coordinate last).   You can use [] for no deformation.
%   options comes from register_phasecorr_initialize. The u_schedule,
%     pixel_spacing, and/or pyramid fields are particularly relevant for
%     this function.
% and
%   sz_unew is the spatial size of the new grid
%   uh is the deformation interpolated to this new grid.
%
% See also: register_phasecorr_initialize.

% Copyright 2010 by Timothy E. Holy

  if iscell(u)
    szu = size(u{1});
  else
    szu = size(u);
  end
  szu(end+1:options.n_dims) = 1;  % ensure that the size vector contains sufficient entries
  szu = szu(1:options.n_dims);  % truncate any n_value component

  %% Calculate the new size
  if any(szu == 0)
    % Input u was empty, start with a global rigid translation
    sz_unew = ones(1,options.n_dims);
    uh = num2cell(zeros(1,options.n_dims));
  else
    % Determine whether this szu is in a lookup table of scheduled
    % refinements
    have_newsz = false;
    if isfield(options,'pyramid')
      if all(szu == 1)
        have_newsz = true;
        sz_unew = options.pyramid(end).sz;
      else
        matched = false;
        for sizeIndex = 2:length(options.pyramid)
          if all(szu == options.pyramid(sizeIndex).sz)
            matched = true;
            break
          end
        end
        if matched
          sz_unew = options.pyramid(sizeIndex-1).sz;
          have_newsz = true;
        else
          error('Using a pyramid, but failed to match the size. Did you already get to the top size?');
        end
      end
    elseif isfield(options,'u_schedule')
      matched = false;
      for sizeIndex = 1:size(options.u_schedule,1)
        if all(szu == options.u_schedule(sizeIndex,:))
          matched = true;
          break
        end
      end
      if matched && (sizeIndex < size(options.u_schedule,1))
        sz_unew = options.u_schedule(sizeIndex+1,:);
        have_newsz = true;
      end
    end
    if ~have_newsz
      if isfield(options,'pixel_spacing')
        % Aim for chunks that have equal physical size along all axes
        sz_pix = options.sz_spatial ./ szu;
        sz_phys = sz_pix .* options.pixel_spacing;
        fac = sz_phys / min(sz_phys);
        fac(fac > sqrt(2)) = 2;
        fac = round(fac);
        if ~any(fac == 2)
          fac(:) = 2;
        end
        sz_unew = fac.*szu - (fac > 1 & szu > 1);
      else
        sz_unew = 2*szu - (szu > 1);
      end
    end
    
    %% Interpolate the old u to the larger size
    if (nargout > 1)
      colons = repmat({':'},1,options.n_dims);
      uh = cell(1,options.n_dims);
      if all(szu == 1)
        for dimIndex = 1:options.n_dims
          if iscell(u)
            thisu = u{dimIndex};
          else
            thisu = u(colons{:},dimIndex);
          end
          uh{dimIndex} = repmat(thisu,sz_unew);
        end
      else
        if isfield(options,'pyramid');
          % Use prolongation (it's much faster)
          for dimIndex = 1:options.n_dims
            if iscell(u)
              thisu = u{dimIndex};
            else
              thisu = u(colons{:},dimIndex);
            end
            uh{dimIndex} = array_prolong(thisu,sz_unew);
          end
        else
          % Use generic interpolation
          % Any dimensions that have size 1 are a problem for interpolation.
          % First, any that will grow, replicate them
          repflag = szu == 1 & sz_unew > 1;
          if any(repflag)
            repsz = ones(1,options.n_dims);
            repsz(repflag) = sz_unew(repflag);
            if iscell(u)
              for dimIndex = 1:options.n_dims
                u{dimIndex} = repmat(u{dimIndex},repsz);
              end
            else
              u = repmat(u,[repsz 1]);
            end
            szu(repflag) = repsz(repflag);
          end
          % Second, any dimensions that will remain of size 1 need to be
          % temporarily discarded
          sz1flag = sz_unew > 1;
          sz1index = find(sz1flag);
          n_dims_i = length(sz1index);
          if (n_dims_i == 1)
            error('Too many unity dimensions. To implement, use interp1')
          end
          x = cell(1,n_dims_i);
          for dimIndex = 1:n_dims_i
            dimIndex2 = sz1index(dimIndex);
            x{dimIndex} = linspace(1,szu(dimIndex2),sz_unew(dimIndex2));
          end
          X = cell(1,n_dims_i);
          [X{:}] = ndgrid(x{:});
          for dimIndex = 1:options.n_dims
            if (sz_unew(dimIndex) > 1)
              if iscell(u)
                uh{dimIndex} = interpn(u{dimIndex},X{:},'linear');
              else
                uh{dimIndex} = interpn(u(colons{:},dimIndex),X{:},'linear');
              end
            else
              uh{dimIndex} = ones(sz_unew,options.class);
            end
          end
        end
      end
    end
  end