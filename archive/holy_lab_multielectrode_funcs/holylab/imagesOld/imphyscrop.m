function ipout = imphyscrop(ipin,cropstruct)
% IMPHYSCROP: crop images in an IMPHYS structure
% Syntax:
%   ipout = imphyscrop(ipin,cropstruct)
% where
%   ipin is the input IMPHYS structure;
%   cropstruct is a structure with the (optional) fields: xrange, yrange,
%     zrange.  Any assigned field gives the [start end] range of
%     retained coordinates along the corresponding dimension; if absent,
%     the entire range is retained.
%
% See also: IMPHYS.
  
% Copyright 2004 by Timothy E. Holy
  ipout = ipin;
  fn = {'xrange','yrange','zrange'};
  for i = 1:length(ipout)
    % Crop, if image is already loaded
    if (isfield(ipout,'image') & ~isempty(ipout(i).image))
      for j = 1:3
        r{j} = ipcrange(ipout(i),cropstruct,fn{j});
      end
      if (ndims(ipout(i).image) == 2)
        ipout(i).image = ipout(i).image(r{1:2});
      else        
        ipout(i).image = ipout(i).image(r{:});
      end
    end
    % Update range info
    for j = 1:3
      if isfield(cropstruct,fn{j})
        ipout(i).(fn{j}) = IntersectIntervals(cropstruct.(fn{j}),...
          ipout(i).(fn{j}));
      end
    end
  end
  
function rng = ipcrange(ip,cropstruct,fn)
  refrng = ip.(fn);
  if isfield(cropstruct,fn)
    rngi = cropstruct.(fn) - refrng(1) + 1;
    if (any(rngi < 1) | any(rngi > diff(refrng)+1))
      error(['Crop on field ' fn ' is out-of-range']);
    end
    rng = rngi(1):rngi(2);
  else
    rng = 1:diff(refrng)+1;
  end
