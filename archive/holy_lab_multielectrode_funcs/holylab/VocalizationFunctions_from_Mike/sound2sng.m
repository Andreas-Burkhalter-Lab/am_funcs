function [clipt] = sound2sng(infilename,sngparms,outfilename)
% SOUND2SNG: calculate the sonogram from the raw waveform
%
% There are usage tips at the end of this help section.
%
% Syntax:
%   sound2sng(infilename,sngparms,outfilename)
%   sound2sng(fid,sngparms,outfilename)
% where
%   infilename is the name of a AI file (.bin extension) OR WAV file (.wav
%     extension)
%   fid is the file identifier of an open AI file (useful if the
%     file needs to opened with a change in endian status, for example.)
%   outfilename is a string containing the name of the output file. If
%     absent, the basename + .sng is chosen (only if a filename rather than
%     a fid is supplied for the input file).
%
% sngparms is a structure with the following fields:
%   threshold:  the threshold (in volt units)
%   nfreq: the number of frequencies in the sonogram
%   freqrange: the range of frequencies (in Hz), [fmin fmax], to keep.
%         Bins which exceed the threshold but do not fall within
%         this frequency range are discarded.
%   plot: if true, gives a semilog plot of the power spectrum.
%         This is useful for inspecting to set the threshold, but
%         not recommended for general use!
%   memmax (optional): the amount of memory to use during processing.
%         With WAV files, you might be able to increase the speed of
%         processing by reading the data in larger chunks. (default: 5Meg)
%   scalemult,scaleoff,voltageMin,voltageMax (optional): if you are reading
%         a WAV file, you can supply these values to convert to microphone
%         voltage. If s is the value stored as a native sample (see
%         WAVREAD), then the voltage is
%                 v = scalemult * s + scaleoff
%         In particular, if the DC-value of your detector is nonzero, you
%         might want to supply scaleoff so you can subtract it.
%         voltageMin and voltageMax store the full-scale range, so you can
%         detect clipping. (defaults: 1, 0, -Inf, Inf, respectively.)
%
% If you're getting started for the first time, the recommended procedure
% is this:
%   1. Find out the noise level of your detector in the units used here:
%        sngparms.plot = true;
%        sngparms.threshold = 0;
%        sngparms.nfreq = 256;
%        sngparms.freqrange = [25000 110000];
%        sound2sng('infile',sngparms,'outfile')
%      Hit control-C to stop the computation. Examine the power spectrum
%      to find the white-noise level (basically, the maximum value
%      attained in the "flat" region of the PSD).
%   2. Now change your parameters to this:
%        sngparms.plot = false;
%        sngparms.threshold = <thresh>;
%      Here <thresh> is whatever you just decided for the white noise
%      level.  The other parameters (nfreq & freqrange) can stay the same
%      (or you can change them if you wish, but there's a chance you
%      might need to tweak thresh again, depending on normalization---I
%      no longer remember. Safest to check.)
%
% See also: SNG2SOUND, WRITESONOGRAM, WAVREAD.
  
% Copyright 2001 & 2008 by Timothy E. Holy <holy@pcg.wustl.edu>

%% Parse the arguments
iswav = false;
if isnumeric(infilename)
  fidSound = infilename;
else
  disp(infilename);
  [fidSound,message] = fopen(infilename,'r');
  if (fidSound == -1)
    error(message);
  end
  [pathstr,name,ext] = fileparts(infilename);
  iswav = strcmp(lower(ext),'.wav');
  if (nargin < 3)
    outfilename = [name '.sng'];
  end
end
disp(['Output file: ',outfilename]);
[fidSng,message] = fopen(outfilename,'w');
if (fidSng == -1)
  error(message);
end
if ~isfield(sngparms,'memmax')
  sngparms.memmax = 5*1024*1024;
end

%% Read the header
if iswav
  sz = wavread(infilename,'size');
  [y,Fs,Nbits] = wavread(infilename,1); % just read 1 sample, so don't overfill memory
  % Fill in header information (fake it where necessary)
  header = struct('nscans',sz(1),'numch',sz(2),'channels',0:sz(2)-1,...
    'scanrate',Fs,'scalemult',1,'scaleoff',0,'voltageMin',-Inf,'voltageMax',Inf,...
    'date','','time','','usrhdr','','flag','WU1');
  header = copyfields(sngparms,{'scalemult','scaleoff','voltageMin','voltageMax'},header);
else
  header = ReadAIHeaderWU1(fidSound);
end

%% Initialize frequency information
zer = zeros(1,sngparms.nfreq);
xold = zer;
columnOffset = 0;

freqRange = [1 sngparms.nfreq];
if isfield(sngparms,'freqrange')
  freqRange = round(sngparms.freqrange/header.scanrate * 2 * (sngparms.nfreq+1));
  freqRange(2) = min([freqRange(2) sngparms.nfreq+1]);
  if (diff(freqRange) <= 0)
    error('Frequency range [fmin fmax] has fmin > fmax');
  end
end
header.nfreq = sngparms.nfreq;

%% Determine timing info and size of blocks to process
header.tacq = header.nscans/header.scanrate;
ncolfactors = factor(floor(header.nscans/(2*sngparms.nfreq)));
nperblock = 1;
% Increase the block size (to boost processing speed), but keep it so
% that there are an integer number of blocks.  We do this by "absorbing"
% factors of the total number of columns.
while (~isempty(ncolfactors) && nperblock*2*sngparms.nfreq*ncolfactors(1) < sngparms.memmax)
  nperblock = nperblock * ncolfactors(1);
  ncolfactors(1) = [];
end
if ~isempty(ncolfactors)
  nblocks = prod(ncolfactors);
else
  nblocks = 1;
end
blockSize = 2 * sngparms.nfreq * nperblock;
header.columnTotal = (blockSize/sngparms.nfreq)*nblocks; % No "2" because half-overlap
header.threshold = sngparms.threshold;
header.nblocks = nblocks;
header.freqMin = freqRange(1);
header.freqMax = freqRange(2);
positionSng = WriteAIHeaderWU1(fidSng,header,struct('AI',1,'Sng',1,'Detect',0));

%% Set up progress bar
if ~isfield(sngparms,'progressbar')
    if(isnumeric(infilename))
        tt_filename=['fid=' num2str(infilename)];
    else
        tt_filename=infilename;
    end
    tName=['File ' tt_filename ': %d %% done'];
    figProgress=waitbar(0,'sound2sng() is busy, please wait...', 'name', sprintf(tName, 0));
elseif isfield(sngparms,'progressbar') && sngparms.progressbar == 1;
    if(isnumeric(infilename))
        tt_filename=['fid=' num2str(infilename)];
    else
        tt_filename=infilename;
    end
    tName=['File ' tt_filename ': %d %% done'];
    figProgress=waitbar(0,'sound2sng() is busy, please wait...', 'name', sprintf(tName, 0));
end

%% Process the data in blocks
fprintf('\n');
dataoffset = 0;
clipt = [];
for i = 1:nblocks
  if (mod(i,floor(nblocks/50)) == 0)
    fprintf('.');
  end
  if iswav
    readbuf = wavread(infilename,[1 blockSize]+dataoffset,'native');
    nread = length(readbuf);
    dataoffset = dataoffset+nread;
  else
    if (isfield(header,'flag') & strcmp(header.flag,'Harvard'))
      [readbuf,nread] = fread(fidSound,[1 blockSize],'int16');
    else
      [readbuf,nread] = fread(fidSound,[1 blockSize],'uint16');
    end
  end
  if (nread < blockSize)
    error('End of file reached prematurely');
    %readbuf = readbuf(1:nread);
  end
  % Calculate actual voltages
  microphoneVoltage = header.scalemult * double(readbuf) + header.scaleoff;
  if any(microphoneVoltage <= header.voltageMin) || any(microphoneVoltage >= header.voltageMax)
    warning('Clipping observed');
    clipt = [clipt;(find(microphoneVoltage <= header.voltageMin |...
        microphoneVoltage >= header.voltageMax)+dataoffset-nread-1)./header.scanrate];
  end
  % Compute the specgram
  [B,xold] = specgramwrap(microphoneVoltage,xold,2*sngparms.nfreq);
  absB = abs(B);
  % Write sonogram data to disk
  [nrows,ncols] = size(B);
  WriteSonogram(fidSng,B,sngparms.threshold,freqRange,columnOffset,absB);
  % Plot to screen?
  if (isfield(sngparms,'plot') && sngparms.plot)
    semilogy(absB);
    drawnow
  end
  columnOffset = columnOffset + ncols;

  % show progress:
  if ~isfield(sngparms,'progressbar')
  if(ishandle(figProgress) && mod(i,floor(nblocks/50)) == 0)
    set(figProgress, 'name', sprintf(tName, round(i/nblocks*100)));
    waitbar(i/nblocks, figProgress);
    % drawnow;
  end
  elseif isfield(sngparms,'progressbar') && sngparms.progressbar == 1
    set(figProgress, 'name', sprintf(tName, round(i/nblocks*100)));
    waitbar(i/nblocks, figProgress);
    % drawnow; 
  end
end
fprintf('\n');
if exist('figProgress','var') && ishandle(figProgress)
    close(figProgress); 
end

%Close the raw data file
status = fclose(fidSound);
if (status < 0)
  error('File did not close');
end
%Save last column to disk
% (needed to re-create the sound from the sonogram)
[B1,xold] = specgramwrap(zer,xold,2*sngparms.nfreq);
WriteSonogram(fidSng,B1,sngparms.threshold,freqRange,columnOffset);
%Close the sparse sonogram file
status = fclose(fidSng);
if (status < 0)
  error('File did not close');
end
