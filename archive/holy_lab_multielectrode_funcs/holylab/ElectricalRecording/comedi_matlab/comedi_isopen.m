function comedi_isopen(devname)
% COMEDI_ISOPEN: Open a data acquisition (comedi) device
% Syntax:
%   ret = comedi_isopen(devicename),
% where
%   devicename is a string giving the device name,
%   ret will be 0 if the device is not open, 1 if it is open.
% This is a MEX file.
%
% See also: COMEDI_OPEN, COMEDI_CLOSE.
