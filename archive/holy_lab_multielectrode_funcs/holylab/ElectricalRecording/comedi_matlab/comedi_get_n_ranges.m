function nranges = comedi_get_n_ranges(devname,subdev,channel)
% COMEDI_GET_N_RANGES: Determine the number of voltage ranges for a channel
% Syntax:
%   nranges = comedi_get_n_ranges(devname,subdev,channel)
% where
%   devname is a string, the name of the device
%   subdev is the subdevice number (see COMEDI_FIND_SUBDEVICE_BY_TYPE)
%   channel is the channel number (starting with 0)
% and
%   nranges is the number of different voltage ranges supported for
%   this channel.
%
% See also: COMEDI_GET_RANGE, COMEDI_FIND_RANGE, COMEDI_FIND_SUBDEVICE_BY_TYPE.
