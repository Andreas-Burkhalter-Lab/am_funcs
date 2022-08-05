function range = comedi_get_range(devname,subdev,channel,rindex)
% COMEDI_GET_RANGE: Determine the voltage range corresponding to range index
% Syntax:
%   [range,unit] = comedi_get_range(devname,subdev,channel,rindex)
% where
%   devname is a string, the name of the device
%   subdev is the subdevice number (see COMEDI_FIND_SUBDEVICE_BY_TYPE)
%   channel is the channel number (starting with 0)
%   rindex is the range number (see COMEDI_GET_N_RANGES)
% and
%   range is a 2-vector of [min max] quoted in units
%   unit is a string describing the units (UNIT_volt, UNIT_mA, or UNIT_none)
%
% See also: COMEDI_GET_N_RANGES, COMEDI_FIND_RANGE, COMEDI_FIND_SUBDEVICE_BY_TYPE.
