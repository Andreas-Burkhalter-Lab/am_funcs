function cmd = cmd_ai_perpetual(subdevice,channels,freq,range)
% CMD_AI_PERPETUAL: make a perpetual AI comedi command
% Syntax:
%   cmd = cmd_ai_perpetual(subdevice,channels,freq,range)
% where
%   subdevice is the number of an analog input device
%   channels is a vector containing the list of channels that you
%     wish to scan, in the proper scan order
%   freq is the scan frequency, in Hz (a scalar)
%   range is the range index, a scalar (used for all channels)
%
% This sets up a command that keeps running until you call
% comedi_cancel_and_flush.  AREF is set to GROUND.  The command starts
% immediately when executed.  You can of course change all these parameters
% by adjusting the fields of the returned cmd structure.
%
% See also: COMEDI_COMMAND, COMEDIBUF_READ, COMEDI_GET_RANGE, CMD_AO_TRIG.
  
  channels = channels(:);
  if (length(freq) > 1)
    error('Frequency must be a scalar');
  end
  if (length(range) > 1)
    error('Range must be a scalar');
  end
  cmd = struct('subdev',subdevice,'flags',0,...
               'start_src','TRIG_NOW','start_arg',0,...
               'scan_begin_src','TRIG_TIMER','scan_begin_arg',1e9/freq,...
               'convert_src','TRIG_TIMER','convert_arg',1e9/freq/length(channels),...
               'scan_end_src','TRIG_COUNT','scan_end_arg',length(channels),...
               'stop_src','TRIG_NONE','stop_arg',0,....
               'chanlist',channels,...
               'rangelist',zeros(size(channels)) + range,...
               'areflist',zeros(size(channels)));

