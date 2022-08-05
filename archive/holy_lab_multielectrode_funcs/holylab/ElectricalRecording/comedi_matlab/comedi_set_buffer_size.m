function comedi_set_buffer_size(devname,subdev,bufsz)
% COMEDI_SET_BUFFER_SIZE: Set the buffer size (for streaming data).
% Syntax:
%   comedi_set_buffer_size(devname,subdev,bufsz)
% where
%   devname is a string giving the device name (must already be open)
%   subdev is the subdevice number (see COMEDI_FIND_SUBDEVICE_BY_TYPE)
%   bufsz is the size of the buffer. (This must be a multiple of the
%     virtual memory page size. See comedilib documentation or just
%     use the value from comedi_get_max_buffer_size.)
%
% See also: COMEDI_GET_MAX_BUFFER_SIZE, COMEDI_FIND_SUBDEVICE_BY_TYPE.
