function rindex = comedi_find_range(devname,subdev,channel,units,range)
% COMEDI_FIND_RANGE: Find the best range index from a range
% Syntax:
%   rindex = comedi_find_range(devname,subdev,channel,units,range)
% where
%   devname is a string, the name of the device
%   subdev is the subdevice number (see COMEDI_FIND_SUBDEVICE_BY_TYPE)
%   channel is the channel number (starting with 0)
%   unit is a string describing the units (UNIT_volt, UNIT_mA, or UNIT_none)
%   range is a 2-vector of [min max] quoted in units
% and
%   rindex is the range number (see COMEDI_GET_N_RANGES)
%
% See also: COMEDI_GET_N_RANGES, COMEDI_GET_RANGE, COMEDI_FIND_SUBDEVICE_BY_TYPE.
