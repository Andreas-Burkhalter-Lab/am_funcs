function ret = comedi_get_max_buffer_size(devname,subdev)
% COMEDI_GET_MAX_BUFFER_SIZE: Determine the maximum allowable buffer size
% Syntax:
%   ret = comedi_get_max_buffer_size(devname,subdev)
% where
%   devname is a string giving the device name (must already be open)
%   subdev is the subdevice number (see COMEDI_FIND_SUBDEVICE_BY_TYPE)
% and
%   ret is the maximum buffer size.
% This max buffer size may be used as an input to
% comedi_set_buffer_size. The maximum buffer size can be adjusted
% only with root permissions; it is often set on bootup using an
% option to comedi_config.
%
% See also: COMEDI_SET_BUFFER_SIZE, COMEDI_FIND_SUBDEVICE_BY_TYPE. 
