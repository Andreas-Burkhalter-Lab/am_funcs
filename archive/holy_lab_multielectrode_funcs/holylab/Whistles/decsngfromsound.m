function sngout = decsngfromsound(infilename,sngparms)
% DECSNGFROMSOUND: compute a decimated sonogram from .bin file
% This works by averaging time bins together. You might be better off
% with the "sparse sonogram" utilities in SOUND2SNG.
%
% Syntax:
%   sngout = decsngfromsound(infilename,sngparms)
% where:
%   infilename is a string giving the .bin file name.  You
%     can also pass a file identifier instead of a file name, if you
%     already have the file open.
%   sngparms is a structure with the following fields:
%     nfreq: the number of frequencies in the sonogram
%     navg:  the number of successive columns to average
% and
%   sngout (optional) is the output decimated sonogram.
%     If absent, the sonogram is plotted in a figure
%
% See also: SOUND2SNG, SNG2SOUND.   

% Copyright 2002 by Timothy E. Holy, <holy@pcg.wustl.edu>
  
  if isnumeric(infilename)
      fidSound = infilename;
  else
    [fidSound,message] = fopen(infilename,'r');
    if (fidSound == -1)
      error(message);
    end
  end

  plotting = 0;
  if (nargout == 0)
    plotting = 1;
  end

  if plotting
    figure
    colormap(1-gray);
  end
  
  header = ReadAIHeader(fidSound);
  tacq = header.nscans/header.scanrate;

  memmax = 5*1024*1024;
  zer = zeros(1,sngparms.nfreq);
  xold = zer;
  columnOffset = 0;
  
  ncolfactors = factor(header.nscans/(2*sngparms.nfreq));
  nperblock = 1;
  while (nperblock*2*sngparms.nfreq*ncolfactors(1) < memmax)
      nperblock = nperblock * ncolfactors(1);
      ncolfactors(1) = [];
  end
  nblocks = prod(ncolfactors);
  blockSize = 2 * sngparms.nfreq * nperblock;

  columnSize = nperblock * nblocks/sngparms.navg;
  sng = zeros(sngparms.nfreq + 1, columnSize);
  nperblock = nperblock/sngparms.navg;
  
  if ~isfield(header,'voltageMax')
      header.voltageMax = 5;
  end
  W = hanning(2*sngparms.nfreq);
  WSS = 2*sngparms.nfreq*sum(power(W,2));
  clim1 = log10((1/4096)*(1/sqrt(2*sngparms.nfreq))*(header.voltageMax)*sqrt(WSS));
  clim2 = log10(0.5*header.voltageMax*sqrt(WSS));

  for nbufs = 0:nblocks-1
    [readbuf,nread] = fread(fidSound,[1 blockSize],'int16');
    if (nread < blockSize)
        error('End of file reached prematurely');
        %readbuf = readbuf(1:nread);    
    end    
    % Calculate actual voltages
    microphoneVoltage = header.scalemult * double(readbuf) + header.scaleoff;
    % Compute the specgram
    [B,xold] = specgramwrap(microphoneVoltage,xold,2*sngparms.nfreq);
    absB = abs(B);
    if (nbufs == 0)
      absB(:,1) = absB(:,2);    % First column contaminated by zero-padding
    end    
    Bo = sngdecimate(absB,2*sngparms.navg);    % 2* because specgramwrap uses half-overlapping windows
    sng(:,((nbufs*nperblock+1):((nbufs+1)*nperblock))) = log10(Bo);
    if plotting
      % Graph sonogram
      imagesc([0 tacq],[0 header.scanrate/2000],sng,[clim1 clim2]);
      axis xy;
      xlabel('Time (s)');
      ylabel('Frequency (kHz)');
      set(gca,'TickDir','out');
      drawnow
    end
  end
  if (nargout > 0)
    sngout = sng;
  end


