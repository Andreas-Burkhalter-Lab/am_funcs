function data = comedi_data_read(devname,subdev,channel,rindex,aref)
% COMEDI_DATA_READ: Read a single point on one channel
% Syntax:
%   data = comedi_data_read(devname,subdev,channel,rindex,aref)
% where
%   devname is a string, the name of the device
%   subdev is the subdevice number (see COMEDI_FIND_SUBDEVICE_BY_TYPE)
%   channel is the channel number (starting with 0)
%   rindex is the range number (see COMEDI_GET_N_RANGES)
%   aref is a string specifying the analog input reference
%     (AREF_GROUND, AREF_COMMON, AREF_DIFF, or AREF_OTHER)
% and
%   data is the value read from the channel (in A/D units).
%
% See also: COMEDI_COMMAND, COMEDI_GET_RANGE, COMEDI_FIND_SUBDEVICE_BY_TYPE.
