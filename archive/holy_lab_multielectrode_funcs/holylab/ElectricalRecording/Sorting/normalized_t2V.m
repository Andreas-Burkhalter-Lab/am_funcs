function [t2V,options] = normalized_t2V(shc,options)
% NORMALIZED_T2V: scale time so it covers span of voltage
%%
% Syntax:
%   t2V = normalized_t2V(shc)
%   [t2V,options] = normalized_t2V(shc,options)
% where
%   shc is a sortheader specialized to a specific channel, see
%     SORTHEAD_IMPORTCHAN;
%   options is a structure which may have the following fields:
%     t2V_tmult: multiplier to apply to the time axis to change its magnitude
%       relative to that of the voltage axes (default: 1);
% and
%   t2V is a scalar, the sought-after conversion factor (t in units of
%     seconds);
%   options is an output options structure with any unsupplied fields
%     filled in with their default values.
%
% See also: ESTIMATE_CLUST_MOVEMENT, SORTHEAD_IMPORTCHAN, AUTOSORT_CALC_T2V.
 
% Copyright 2005 by Timothy E. Holy

  if (nargin < 2)
    options = struct;
  end
  if ~isfield(options,'t2V_tmult')
    options.t2V_tmult = 1;
  end

  vmax = -Inf;
  vmin = Inf;
  for i = 1:length(shc)
    peakheight = double(shc(i).peakHeight)*shc(i).scalemult + shc(i).scaleoff;
    vmax = max(vmax,max(peakheight));
    vmin = min(vmin,min(peakheight));
  end
  t = sortheader_absolute_sniptime(shc,'all');
  if iscell(t)
    t = cat(2,t{:});
  end
  t2V = options.t2V_tmult * (vmax-vmin)/(t(end)-t(1));
  