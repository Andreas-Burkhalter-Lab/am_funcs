function comedi_dio_config(devname,subdevice,channel,direction)
% COMEDI_DIO_CONFIG: configure digital I/O channels for input or output
% Syntax:
%    comedi_dio_config(devname,subdevice,channel,direction)
% where
%   devname is a string, the name of the device
%   subdevice is an integer, the number of the subdevice (see
%     COMEDI_FIND_SUBDEVICE_BY_NAME)
%   channel is the digital channel number (starting at 0)
%   direction is a string which may take the values COMEDI_INPUT or
%     COMEDI_OUTPUT. 
%
% See also: COMEDI_OPEN, COMEDI_FIND_SUBDEVICE_BY_TYPE.
