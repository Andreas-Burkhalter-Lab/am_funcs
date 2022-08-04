function htb = htbInit;
% htbInit - Initialize HTB structure
% Returns an empty HTB structure.
%
% SYNOPSIS
%       htb = htbInit;
%
% IN
%       nothing
%
% OUT
%       An empty HTB structure.
%       All numeric elements are 0 except htb.fid which is -1.
%       All string elements are ''.

htb.date                    = '';
htb.ldate                   = 0;
htb.cfg_file                = '';
htb.pro_file                = '';
htb.speed                   = 0;
htb.offset                  = 0;
htb.period                  = 0;
htb.extension               = 0;
htb.skip                    = 0;
htb.first_channel           = 0;
htb.nchannels               = 0;
htb.sweep_limit             = 0;
htb.cancel_override         = 0;
htb.func                    = 0;
htb.tag                     = 0;
htb.npages                  = 0;
htb.nsamples                = 0;
htb.samples_per_page        = 0;
htb.sweep                   = 0;
htb.next_page               = 0;
htb.next_off                = 0;
htb.title                   = '';
htb.speed_units             = 0;

htb.fileoffset = 0;         % file offset for later use
htb.fid = -1;               % file id for later use

return;
