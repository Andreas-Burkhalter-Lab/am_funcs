function [status,dataout] = SonySerialGate(command,datain)
% SonySerialGate: drive the camcorder from MATLAB
% [status,dataout] = SonySerialGate(command,datain)
%
% The command input drives the camcorder according to the following table:
%                0        Stop
%                1        Pause
%                2        Play
%                3        Record
%                4        ReadTimeCode (requires 2 outputs: timecode is returned in dataout)
%                5        Cue to timecode (requires 2 inputs, with the timecode in datain)
%                6        CueReady (is the cueing operation complete?)
%                7        ScanForward (FF while viewing on screen)
%                8        ScanBackward
%                9        SlowForward (play at 1/10 speed)
%                10        SlowBackward (play backwards at 1/10 speed)
%                11        Rewind
%                12        Fast forward
%
% This is written as a MEX file.
