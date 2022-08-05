function wout = wiener_resize(win,len)
% WIENER_RESIZE: increase the size of a Wiener filter.
% This basically handles the wrapping issues.
% Syntax:
%   wout = wiener_resize(win,len)
% where
%   win is the input filter (e.g., from WIENER).
%   len is the desired size
% and
%   wout is the output resized vector.
%
% See also: WIENER.

% Copyright 2007 by Timothy E. Holy

  [clen,n_channels,n_templates] = size(win);
  wout = zeros([len,n_channels,n_templates],class(win));
  partdiv = ceil(clen/2);
  wout(1:partdiv,:,:) = win(1:partdiv,:,:);
  wout(end-(clen-partdiv)+1:end,:,:) = win(partdiv+1:end,:,:);
  