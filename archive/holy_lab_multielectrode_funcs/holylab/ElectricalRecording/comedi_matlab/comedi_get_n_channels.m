function nchannels = comedi_get_n_channels(devname,subdev)
% COMEDI_GET_N_CHANNELS: Determine the number of channels in subdevice
% Syntax:
%   nchannels = comedi_get_n_channels(devname,subdev)
% where
%   devname is a string, the name of the device
%   subdev is the subdevice number (see COMEDI_FIND_SUBDEVICE_BY_TYPE)
% and
%   nchannels is the number of different voltage channels supported for
%   this channel.
%
% See also: COMEDI_GET_RANGE, COMEDI_FIND_RANGE, COMEDI_FIND_SUBDEVICE_BY_TYPE.
