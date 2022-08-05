function success=comedi_internal_trigger(devname, subdevice, trignum)
% comedi_internal_trigger: initiate a waiting command
%
% A command can be prepared with start_src = 'TRIG_INT'; when executed,
% this command waits until it receives a trigger.  This function sends
% such a trigger.
%
% Syntax:
%   success=comedi_internal_trigger(devname, subdevice, trignum)
% where
%   devname is a string, the name of the device
%   subdevice is the subdevice number,
%   trignum is the trigger number.
% and
%   success is <0 if this command fails.
%
% SEE ALSO: comedi_command, comedibuf_write, comedi_find_subdevice_by_type.
