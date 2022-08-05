function comedi_cancel_and_flush(devname,subdev)
% COMEDI_CANCEL_AND_FLUSH: Terminate a board-autonomous command
% This is useful when you want to abort, or if you set up your
% command with no stop condition.
% Syntax:
%   comedi_cancel_and_flush(devname,subdev)
% where
%   devname is a string, the name of the device
%   subdev is the subdevice number (see COMEDI_FIND_SUBDEVICE_BY_TYPE)
%
% For a subdevice of type AI, the read buffer is also read out (and dumped)
% until it empties.  Hopefully, this prepares the subdevice for
% further action.
%
% See also: COMEDI_COMMAND, COMEDI_FIND_SUBDEVICE_BY_TYPE.
