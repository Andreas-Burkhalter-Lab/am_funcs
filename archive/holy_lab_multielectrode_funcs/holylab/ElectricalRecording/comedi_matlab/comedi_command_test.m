function [msg,cmdo] = comedi_command_test(devname,command)
% COMEDI_COMMAND_TEST: Check a board-autonomous command
% Syntax:
%   [msg,cmdo] = comedi_command_test(devname,command)
% where
%   devname is a string, the name of the device
%   command is a MATLAB structure which specifies the command.
%     See MAKE_COMEDI_CMD for an example, and COMEDI_COMMAND for a
%     description.
% and
%   msg is a string which tells you the outcome of the test.
%     'success' indicates that everything is ready to go.
%   cmdo is a modified version of your original command, tweaked to
%     fit the supported values for your board.
%
% Frequently, a user-supplied command will need to be passed
% through COMEDI_COMMAND_TEST several times before all the tweaking
% is done.  You can embed the call to COMEDI_COMMAND_TEST within a
% while loop, checking for msg == 'success'.
%
% See also: MAKE_COMEDI_CMD, COMEDI_COMMAND.
