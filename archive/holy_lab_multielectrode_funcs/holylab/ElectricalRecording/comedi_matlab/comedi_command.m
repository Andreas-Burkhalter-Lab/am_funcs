function ret = comedi_command(devname,command)
% COMEDI_COMMAND: Initiate a board-autonomous command
% Syntax:
%   ret = comedi_command(devname,command)
% where
%   devname is a string, the name of the device
%   command is a MATLAB structure which specifies the command.
%     See MAKE_COMEDI_CMD for an example.
%     The following is largely copied from the "cmd" demo of
%     Comedilib, written by David A. Schleef.
%     Each event requires a trigger, which is specified by a source
%     and an argument.  For example, to specify an external digital
%     line 3 as a source, you would use src=TRIG_EXT and arg=3.
%
%     This command structure has the following fields:
%       subdev:         the subdevice number (see COMEDI_FIND_SUBDEVICE_BY_TYPE)
%       flags:          some flags (see COMEDI_MISSING)
%       start_src:      Controls the start of acquisition. A string
%                       which may take the following values:
%          TRIG_NOW:     The start_src event occurs start_arg nanoseconds
%                        after comedi_command() is called.  Currently,
%                        only start_arg=0 is supported.
%          TRIG_FOLLOW:  (For an output device.)  The start_src event occurs
%                        when data is written to the buffer.
%          TRIG_EXT:     start event occurs when an external trigger
%                        signal occurs, e.g., a rising edge of a digital
%                        line.  start_arg chooses the particular digital
%                        line.
%          TRIG_INT:     start event occurs on a Comedi internal signal,
%                        which is typically caused by an INSN_TRIG
%                        instruction.
%       start_arg:      A scalar, usually 0 (see start_src).
%       scan_begin_src: Controls the timing of the beginning of
%                       each scan.  This is a string with the
%                       following values:
%          TRIG_TIMER:   scan_begin events occur periodically.
%                        The time between scan_begin events is
%                        scan_begin_arg nanoseconds.
%          TRIG_EXT:     scan_begin events occur when an external trigger
%                        signal occurs, e.g., a rising edge of a digital
%                        line.  scan_begin_arg chooses the particular digital
%                        line.
%          TRIG_FOLLOW:  scan_begin events occur immediately after a scan_end
%                        event occurs.
%       scan_begin_arg: See scan_begin_src.  For regular scan
%                       events occuring at a frequency of freq (in
%                       seconds), you would set this to 1e9/freq,
%                       with scan_begin_src set to TRIG_TIMER.  The
%                       scan_begin_arg that we use here may not be
%                       supported exactly by the device, but it
%                       will be adjusted to the nearest supported
%                       value by comedi_command_test.
%
%       convert_src:    Controls the timing between samples in each
%                       scan.  This is a string which may assume
%                       the following values:
%          TRIG_TIMER:   Conversion events occur periodically.
%                        The time between convert events is
%                        convert_arg nanoseconds.
%          TRIG_EXT:     Conversion events occur when an external trigger
%                        signal occurs, e.g., a rising edge of a digital
%                        line.  convert_arg chooses the particular digital
%                        line.
%          TRIG_NOW:     All conversion events in a scan occur simultaneously.
%       convert_arg:    See convert_src.  You can start with 1
%                       (nanosecond), and then comedi_command_test
%                       will adjust this value to a board-supported value.
%
%       scan_end_src:   Controls the end of each scan.  The end of
%                       each scan is almost always specified using
%                       TRIG_COUNT, with the argument being the
%                       same as the number of channels in the
%                       chanlist.  You could probably find a device
%                       that allows something else, but it would be
%                       strange.
%       scan_end_arg:   Usually, the number of channels in the
%                       channel list.
%
%       stop_src:       Controls the end of acquisition.  This is a
%                       string which may assume the following values:
%          TRIG_COUNT:   stop acquisition after stop_arg scans.
%          TRIG_NONE:    continuous acquisition, until stopped using
%                        comedi_cancel_and_flush.
%       stop_arg:       Set to the number of scans, when using
%                       stop_src = 'TRIG_COUNT'.
%
%     Information about the individual channels is supplied in
%     either of two formats:
%     Format #1:
%       chanlist:       A vector of the channel numbers, in the
%                       order in which they are to be scanned
%       rangelist:      The range index for each channel (see
%                       COMEDI_GET_RANGE and associated
%                       functions).  For multiplexed boards, you
%                       may want to use the same range for all
%                       channels, as resetting the gain on the
%                       internal amplifier can take a long time to
%                       settle.
%       areflist:       A vector of integers giving the analog
%                       reference for each channel.  The
%                       permissible values are 0 = AREF_GROUND,
%                       1 = AREF_COMMON, 2 = AREF_DIFF, 3 = AREF_OTHER.
%     Format #2:
%       packedchanlist: A vector, where each entry has been
%                       "packed" by comedi to encode the channel
%                       number, range, and analog reference.  In
%                       general, you probably won't manipulate
%                       this directly; but this form is returned by
%                       functions such as comedi_command_test.
%     You can't use both formats in the same command structure.
%
% It's strongly recommended that you first pass your command
% through comedi_command_test, to see if it's ready to go.  This
% procedure will also tweak values to ones supported by the board.
%
% See also: MAKE_COMEDI_CMD, COMEDI_COMMAND_TEST, COMEDIBUF_READ,
%   COMEDI_CANCEL_AND_FLUSH, COMEDI_DATA_READ, COMEDI_GET_RANGE,
%   COMEDI_FIND_SUBDEVICE_BY_TYPE, COMEDI_MISSING.
