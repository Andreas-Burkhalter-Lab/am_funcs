function [sng,f,t] = sparsesng(y,Fs,NFFT,thresh,options)
% SPARSESNG: sparse representation of the sonogram
% sng = sparsesng(y,Fs,NFFT,thresh)
% OR
% sng = sparsesng(filename,trange,NFFT,thresh)
% OR
% sng = sparsesng(fid,trange,NFFT,thresh)
%  
% where
%    y is a vector signal, and
%    Fs is the sampling frequency
% OR
%    filename is the name of the .bin file (AI type) OR
%    fid is a file identifier of such a file,     and
%    trange is the [starttime endtime], in seconds, to use (if empty: all)
% Also:
%    NFFT (power of 2) is the size of the fourier window
%    thresh is the threshold for retaining values
% and
%    sng is a sparse matrix of threshholded fft components (NFFT/2+1 by ?)
%    (sng only contains the non-redundant half of the fft components, as
%    signals are real.) By default, each block half-overlaps with
%    its neighbor. The data are windowed before taking the fft.
%
% [sng,f,t] = sparsesng(...) returns frequency and time vectors for
%    scaling the "axes" of sng. f is in units of kHz, t is in seconds.
%
% sng = sparsesng(...,options) allows one to set a variety of other
%   parameters through the structure 'options', with fields:
%     band: the 2-vector giving the frequency range (in kHz) to
%       retain.  If not specified, the default is to use [35 100].
%     numoverlap: specifies the number of samples by which
%       successive blocks overlap. The default is NFFT/2.
%     totempfile: dump output to a temporary file, instead of
%       holding in memory. Specify the name of the file as a
%       string. This is currently WISH LIST status.
%
% See also: SNG2SOUND

% Copyright 05-03-01 by Tim Holy
% Changelog 11-01-01 by Tim Holy:
%   Window each block with a hanning window (required changes in
%     SNG2SOUND)
%   Convert to volts, not A/D units
%   Allow fid input (so can set endian flag ahead of time)
% Changelog 11-07-01 by Tim Holy:
%   Make last argument an options structure
%   Allow arbitrary overlaps (needed for interfacing with
%     phasevocoder)
%   Convert band to a 2-vector only (used to be a 4-vector,
%     allowing a control band, but this did not seem useful).
%   Made the time axis (t output) more accurate---it gives
%     the start time of each block
  if (nargin < 5)
    options = [];
  end
  if ~isfield(options,'band')
    options.band = [35 100];
  end
  if ~isfield(options,'numoverlap')
    options.numoverlap = NFFT/2;
  end
  fromfile = 0;
  if (ischar(y))
    fromfile = 1;
  end
  if (isnumeric(y) & length(y) == 1)
    fromfile = 1;   % assume an fid
    fseek(y,0,'bof'); % rewind the file
  end
  tovolts = 1;
  if (fromfile)
    filename = y;
    trange = Fs;
    if ischar(filename)
      % Open file and read header
      [fid,message] = fopen(filename,'r');
      if (fid < 1)
        error(message)
      end
    else
      fid = y;
    end
    h = ReadAIHeader(fid);
    status = fseek(fid,h.headersize,'bof');        % Move to the beginning of data
    if (status)
      error(ferror(fid))
    end
    % Do error checking here
    if (h.numch ~= 1)
      error('Must have exactly 1 recorded channel');
    end
    % Get to start time
    if (isempty(trange))
      trange = [0 h.nscans-1];
    else
      trange = round(trange*h.scanrate);
    end
    status = fseek(fid,trange(1)*2,'cof');    % *2 because int16s
    % Set up generic data
    Fs = h.scanrate;
    npts = diff(trange)+1;
    tovolts = h.scalemult;
  else
    npts = length(y);
  end
  l2 = log2(NFFT);
  if (floor(l2) ~= l2)
    error('NFFT must be a power of 2');
  end
  if (NFFT < options.numoverlap)
    error('Cannot overlap by more than NFFT');
  end
  thewindow = hanning(NFFT)';
  % Compute the mapping between frequency indices and frequencies
  nfreqs = NFFT/2+1;
  nnew = NFFT - options.numoverlap;
  % blocksize = NFFT/2;
  f = linspace(0,Fs/2000,nfreqs);        % bands measured in kHz
  bandsig = find(f >= options.band(1) & f <= options.band(2));
  % Compute the size of the output
  nblocks = floor((npts-options.numoverlap)/nnew);
  indx1 = cell(1,nblocks);
  sng1 = cell(1,nblocks);
  % Load the initial data
  s = zeros(1,NFFT);
  if (fromfile)
    s(nnew+1:end) = tovolts * ReadBinaryData(fid,1,[0 options.numoverlap-1]);
  else
    s(nnew+1:end) = y(1:options.numoverlap);
  end
  % Loop through the data in blocks
  for k = 1:nblocks
    % First, transfer the old data to the overlap region
    s(1:options.numoverlap) = s(nnew+1:end);
    % Now get new data
    if (fromfile)
      s(options.numoverlap+1:end) = tovolts * ReadBinaryData(fid,1,[0 nnew-1]);
    else
      s(options.numoverlap+1:end) = y(options.numoverlap + ...
                                      ((k-1)*nnew+1:k*nnew));
    end
    % Do fourier transform
    st = fft(thewindow .* s,NFFT);
    indx = find(abs(st(bandsig)) > thresh);
    indx1{k} = bandsig(indx);
    sng1{k} = st(bandsig(indx));
    %fprintf('%d out of %d\n',k,nblocks);
  end
  % Now turn this cell-array organization into a sparse matrix
  lennz = zeros(1,nblocks);
  for k = 1:nblocks
    lennz(k) = length(indx1{k});
  end
  indx2 = zeros(1,sum(lennz));
  start = 1;
  for k = 1:nblocks
    indx2(start:start+lennz(k)-1) = k;
    start = start+lennz(k);
  end
  sng = sparse(cat(2,indx1{:}),indx2,cat(2,sng1{:}),nfreqs,nblocks);
  % Temporal coordinates
  t = linspace(0,npts-options.numoverlap,nblocks)/Fs;
