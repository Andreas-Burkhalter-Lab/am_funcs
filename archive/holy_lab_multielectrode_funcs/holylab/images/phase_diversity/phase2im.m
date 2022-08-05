function [im,cb] = phase2im(phi,options)
% PHASE2IM: convert a phase function to an image for plotting
% Syntax:
%   im = phase2im(phase)
%   im = phase2im(phase,options)
% where
%   phase is a matrix indicating the phase as a function of position.
%     You can use NaN's to mark "pixels" that you want to plot in the
%     background color.
%   options is a structure which may have the following fields:
%     cm (default jet): colormap to use
%     bg (default [1 1 1]): background color to use
%     scale: if present, scales the phase by this value (otherwise the
%       maximum absolute value is used)
%     wrap (default false): if true, calls filter_wrap to convert from a
%       periodic "corners"-based phase to a "center"-based phase.
%     truncate (default true): if true, excludes rows/columns that are all
%       NaN.
% and
%   im is an output RGB image suitable for plotting with image()
  
% Copyright 2008 by Timothy E. Holy
  
  if (nargin < 2)
    options = struct;
  end
  options = default(options,'cm',jet,'bg',[1 1 1],'wrap',false,'truncate',true);
  
  if ~isfield(options,'scale')
    % Scale the phase to maximum value
    phimax = max(abs(phi(:)));
    if (phimax == 0)
      phimax = 1;
    end
  else
    phimax = options.scale;
  end
  
  if options.wrap
    phi = filter_wrap(phi,size(phi));
  end
  
  if options.truncate
    nancheck = ~isnan(phi);
    rowkeep = any(nancheck,1);
    colkeep = any(nancheck,2);
    phi = phi(rowkeep,colkeep);
  end
  
  % Calculate the phase color lookup
  cm = options.cm;
  ncm = size(cm,1);
  colindx = round((phi+phimax)/(2*phimax) * (ncm - 1) + 1);
  % Make the NaN pixels the color of the background
  colindx(isnan(phi)) = ncm+1;
  cm(ncm+1,:) = options.bg;
  im = reshape(cm(colindx(:),:),[size(phi) 3]);
end
