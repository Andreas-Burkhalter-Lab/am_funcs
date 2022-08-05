function sop = snipoptions(sopin)
% SNIPOPTIONS: options to control behavior during snippet cutting
% sop = snipoptions(sopin),
% where sop & sopin are structures with the following fields:
%   condfilta,condfiltb (default 1): conditioning filters for the raw
%      waveform (see FILTER). A value of 1 for these results in no filtering.
%   detfilt (default empty): matrix of detection filters, each column
%      corresponding to a channel
%   polarity (default 1): set to 1 for positive-going peaks, -1 for
%      negative-going, and 0 for both
%   close (default 30): number of samples between peaks before they are
%          tested for adjacency
%      troughdepth (default 0.5): factor by which waveform between two adjacent
%          peaks must dip for them to be considered different events
%          (actually - lowest dip must pass below troughdepth*lowerpeak)
%      peaktrough (default 30): minimum number of samples between
%          opposite-going triggers (used only if polarity = 0; prevents triggering
%          on both the peak and the trough of the waveform)
%      reboundS (default 0; 100 is good start point): maximum # samples between 
%          peaks that will trigger them to be screened for actually representing 
%          the rebound voltage of an earlier spike (aka, equivelent of "close"
%          but for rebound screen; if set to 0, rebound screen turned off)
%      reboundC (default 4000): minimum curvature that will trigger a spike to
%          be screened for actually representing the rebound voltage of an
%          earlier spike (if set to 0, rebound screen is turned off)
%      reboundO (default 1): selects which combination of outputs will be 
%          displayed and/or saved from the rebound screen. Saved output 
%          is called basefileName.merec.mat, and consists of a matrix named
%          reboundStats (where the first row contains the curvature values
%          of all spikes that appear to be the third of a rebound-like
%          triad, the second contains the time intervals between the 
%          first and third spike in the triad, the third contains 
%          a 0 if the curvature value was less than reboundC and a 1 otherwise,
%          and the fourth contains the scan number of the third member of the triad)
%          as well as soptions; displayed output consists of two histograms,
%          one of the curvature values and the other of the time intervals, 
%          for all rebound-like triads detected (if multiple files are processed at 
%          once, a separate figure will be made for each file).
%               0 = no output saved or displayed
%               1 = output saved but not displayed
%               2 = output displayed but not saved
%               3 = output saved and displayed
%   blocksize (default 50000): the number of scans to handle at a time
%   interptimes (default 0): if true, causes the output snippet times to
%      be interpolated to a fractional scan number
%   interpsnips (default 0): if true, causes the snippet waveforms to
%      be interpolated to their values at fractional scan numbers
%   tempfilename (default sniptemp.tmp): the name of the temporary file
%      (used only if tofile is true)
%   tofile (default 0): a flag (0 or 1) which, if true, results in the
%      times and snippets being written to a file
%   outfilename (default non-existent): the name of the output file
%      (used only if tofile is true)
%      
%
% The difference between conditioning filters and detection filters
% is that the output snippet waveforms will have been passed through
% the conditioning but not the detection filters; detection filters are
% used only for determining threshold crossing.
%
% The following fields may be set
% in sopin, and SNIPOPTIONS calculates the appropriate conditioning filters:
%   Fs (no default value): must be provided to calculate the conditioning
%      filters, if the fields below are supplied
%   Hz60: a flag (0 or 1) which, if true, sets a notch filter
%      for 60Hz noise
%   bandpass: a 2-vector giving the frequency band
%      (in fractions of the Nyquist frequency) to pass in conditioning
%      the signal
% These are used in addition to whatever is already supplied for conditioning
% filters. On output, these are set to the default values, so that
% successive calls to SNIPOPTIONS do not yield additional filtering. The
% input values are appended to the fields histHz60 and histbandpass,
% respectively.
% If Fs is not supplied, then setting the fields Hz60 and bandpass
% results in an error.
%
% See also: SNIPPETFILE, CUTSNIPPETS, FILTER

% History:
%   2004-08-09: 1. Rebound criteria introduced
%               (RCH)

if (nargin > 0)
  sop = sopin;
else
  sop = [];
end
% Set defaults
if (~isfield(sop,'condfilta') | isempty(sop.condfilta))
  sop.condfilta = 1;
end
if (~isfield(sop,'condfiltb') | isempty(sop.condfiltb))
  sop.condfiltb = 1;
end
if (~isfield(sop,'detfilt') | isempty(sop.detfilt))
  sop.detfilt = [];
end
if (~isfield(sop,'polarity') | isempty(sop.polarity))
  sop.polarity = 1;
end
if (~isfield(sop,'close') | isempty(sop.close))
  sop.close = 30;
end
if (~isfield(sop,'troughdepth') | isempty(sop.troughdepth))
  sop.troughdepth = 0.5;
end
if (~isfield(sop,'peaktrough') | isempty(sop.peaktrough))
  sop.peaktrough = 30;
end
if (~isfield(sop,'reboundS') | isempty(sop.reboundS))
    sop.reboundS = 0;
end
if (~isfield(sop,'reboundC') | isempty(sop.reboundC))
    sop.reboundC = 4000;
end
if (~isfield(sop,'reboundO') | isempty(sop.reboundC))
    sop.reboundO = 1;
end
if (~isfield(sop,'interptimes') | isempty(sop.interptimes))
  sop.interptimes = 0;
end
if (~isfield(sop,'interpsnips') | isempty(sop.interpsnips))
  sop.interpsnips = 0;
end
if (~isfield(sop,'blocksize') | isempty(sop.blocksize))
  sop.blocksize = 50000;
end
if (~isfield(sop,'tofile') | isempty(sop.tofile))
  sop.tofile = 0;
end
if (~isfield(sop,'tempfilename') | isempty(sop.tempfilename))
  if isunix
    sop.tempfilename = '/tmp/sniptemp.tmp';
  else
    sop.tempfilename = 'sniptemp.tmp';
  end
end
if (sop.tofile & (~isfield(sop,'outfilename') | isempty(sop.outfilename)))
  error('Must supply outfilename if saving to file');
end

% Now provide for simple conditioning filtering
if (isfield(sop,'Hz60') & isfield(sop,'Fs'))
  Hz60 = sop.Hz60;
  if (isfield(sop,'histHz60'))
    sop.histHz60(end+1) = sop.Hz60;
  else
    sop.histHz60 = sop.Hz60;
  end
  sop = rmfield(sop,'Hz60');
  if (Hz60)
    %[b,a] = ellip(3,0.5,10,[55 65]/(sop.Fs/2),'stop');
    disp(['Warning: potential stability issues, see help for' ...
          ' filterint16.']);
    disp('Consider using a bandpass instead.');
    [b,a] = butter(2,[50 70]/(sop.Fs/2),'stop');
    sop.condfilta = conv(sop.condfilta,a);
    sop.condfiltb = conv(sop.condfiltb,b);
  end
elseif (isfield(sop,'Hz60') & ~isfield(sop,'Fs'))
  error('Fs must be supplied if Hz60 is set');
end
if (isfield(sop,'bandpass') & ~isempty(sop.bandpass) &  isfield(sop,'Fs'))
  if (length(sop.bandpass) ~= 2)
    error('bandpass must be a 2-vector');
  end
  if (isfield(sop,'histbandpass'))
    sop.histbandpass(end+1,:) = sop.bandpass(:)';
  else
    sop.histbandpass = sop.bandpass(:)';
  end
  [b,a] = butter(2,sop.bandpass/(sop.Fs/2));
  sop = rmfield(sop,'bandpass');
  sop.condfilta = conv(sop.condfilta,a);
  sop.condfiltb = conv(sop.condfiltb,b);
elseif (isfield(sop,'bandpass') & ~isfield(sop,'Fs'))
  error('Fs must be supplied if bandpass is set');
end
