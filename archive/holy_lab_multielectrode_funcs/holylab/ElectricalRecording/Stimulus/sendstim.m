function sendstim(sequence,outfilename)
% SENDSTIM: start a process to run the stimulus sequence
%
% Syntax:
%   sendstim(sequence)
%   sendstim(sequence,outfilename)
% where
%   sequence is a 2-by-ntransitions matrix, with the top row giving
%     the valve numbers and the bottom row the times (monotonic
%     increasing).
%   outfilename is a string giving the name of an optional output file
%
% Examples: if sequence is
%    0    1   0    2
%    0    1.5 1.6  3
% then
%    valve 0 will be turned on right away
%    valve 1 will turn on at t = 1.5sec
%    valve 0 will turn on at t = 1.6sec
%    valve 2 will turn on at t = 3sec, and will be left on.
%
% This sequence runs without blocking MATLAB.
  
  if (nargin < 2)
    outfilename = '/dev/null';
  end
  stimpath = '/usr/lab/matlabfunc/ElectricalRecording/Stimulus/';
  % Make sure there are no stimulus sequences running already
  [status,result] = unix([stimpath 'killsendstim']);
  if (nargin > 0)
    if (size(sequence,1) ~= 2)
      error('sequence input must be a 2-by-n matrix');
    end
    % Prepare the command string
    inptstring = [sprintf('%d\n',size(sequence,2)), ...
                  sprintf('%d %g\n',sequence)];
    cmdstr = ['echo "' inptstring '" | ' stimpath 'sendstim > ' outfilename ...
              ' &'];
    unix(cmdstr);
  end
                
                
