% Only part of the functionality of comedilib has been wrapped for
% MATLAB.  This file documents some of these shortcomings and how
% get around some of them.
%
% The flags field of the command structure must be set with
% integers; none of the nice ORing possible with C is implemented.
% See cmd.c in the comedilib/demo directory for help.
%
% The insn support is very limited: only WAIT and INTTRIG are
% supported.  However, it shouldn't take much to implement other
% instructions (I've just never needed them).  However, with recent
% changes to triggering, I recommend using comedi_internal_trig
% instead of the general INSN framework.
