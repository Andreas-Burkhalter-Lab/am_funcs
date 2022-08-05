function spikes = shape_spikes_to_array(spikes_in,channel_list, ...
					array_channels)
% SHAPE_SPIKES_ON_ARRAY: reshape spike waveform to array geometry
% Syntax:
%   spikes_out = shape_spikes_to_array(spikes_in,channel_list,array_channels)
% where
%   spikes_in is a T-by-n_channels-by-n_spikes array
%      (T = duration of spike waveform,
%       n_channels = # of electrodes
%       n_spikes is the number of overlapping spikes to show)
%   channel_list is a vector of channel #s, in the order used for the 2nd
%     coordinate of "spikes_in"
%   array_channels is a matrix of size ny-by-nx, giving the channel order
%     (see get_hda_holylab).
% and
%   spikes_out is an array of the shape expected by plot_spikes_on_array.
%
% See also: GET_HDA_HOLYLAB, PLOT_SPIKES_ON_ARRAY.

% Copyright 2007 by Timothy E. Holy

  [ny,nx] = size(array_channels);
  [T,n_channels,n_spikes] = size(spikes_in);
  if (n_channels ~= nx*ny)
    error(['The number of array channels is not equal to the number of' ...
	   ' specified channels']);
  end
  try
    p = findainb(array_channels(:),channel_list(:));
  catch
    error('The channel lists do not match');
  end
  spikes = spikes_in(:,p,:);
  spikes = reshape(spikes,[T,ny,nx,n_spikes]);