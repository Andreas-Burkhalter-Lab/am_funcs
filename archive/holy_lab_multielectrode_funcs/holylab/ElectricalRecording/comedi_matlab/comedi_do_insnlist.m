function success = comedi_do_insnlist(devname,insns)
% COMEDI_DO_INSNLIST: execute an instruction list
% Instructions provide one way to execute synchronous tasks in
% comedi.  For devices which support streaming commands,
% instruction lists provide a way to trigger the execution of
% commands that have been set up with start_src = 'TRIG_INT'.
%
% This function provides limited support for comedi's insn
% interface.  Currently, only WAIT and INTTRIG are supported,
% although some of the infrastructure for other instructions is in
% place.
%
% Syntax:
%   success = comedi_do_insnlist(devname,insns)
% where
%   devname is a string, the name of the device
%   insn is a structure array.  Multiple instructions can be sent
%   as different elements of the array.  Here are the different
%   commands available and the fields that must be set to implement
%   them:
%   WAIT:
%     .insn = 'INSN_WAIT'
%     .n = 1
%     .data = # of nanoseconds to wait (up to 10000 is OK, gives me
%        a segfault)
%   INTTRIG:  (send a trigger to a subdevice)
%     .insn = 'INSN_INTTRIG'
%     .subdev = subdevice #
%
% See also: COMEDI_COMMAND, COMEDI_OPEN.
