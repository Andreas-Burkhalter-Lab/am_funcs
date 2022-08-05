function spsngplot(sng,varargin)
% SPSNGPLOT: display a sparse sonogram (uses decimation)
% Usage:
%   spsngplot(filename)  (file must be a sparse sng file)
%   spsngplot(sng)
%   spsngplot(sng,f,t)   (to specify frequency and time axes)
%   spsngplot(...,options)
% where
%   filename is the name of a .sng file
%   sng is the sparse sonogram
%   f is the vector of frequencies
%   t is the vector of times
%   options is an optional structure with the following fields possible:
%     WithAmplitude: if true, plots a second axis containing the
%       amplitude vs. time
%     sliderwindow (default true): if true, creates a sliderwindow.

% Copyright 2002 by Timothy E. Holy <holy@pcg.wustl.edu>

  
  % Argument parsing and defaults
  if ischar(sng)
    filename = sng;
    [sng,header] = ReadSonogram(filename);
    f = linspace(0,header.scanrate/2000,header.nfreq);
    % f sets the y axis.
    % for the sonogram, this is frequency in kHz.
    % linspace sets up linearly spaced on (a,b,n), where a and b
    % are the limits. in this case that's 0 to 125 kHz. header.scanrate
    % is the sampling rate. That was 250Khz. This allows resolution of
    % frequencies up to 125 kHz (Nyquist), thus 250000/2000=125.
    % nfreq is the number of frequencies for the FFT as specified in 
    % sound2sng under the parameter nfreq. The default is 256. 
    
    t = linspace(0,header.nscans/header.scanrate,header.columnTotal);
    
    % So, similarly, this sets up the time axis (x axis) for the sonogram
    

    nextarg = 1;
  else
    nextarg = 1;
    if (length(varargin) >= nextarg & ~isstruct(varargin{nextarg}))
      f = varargin{nextarg};
      nextarg = nextarg + 1;
    else
      f = 1:size(sng,1);
    end
    if (length(varargin) >= nextarg & ~isstruct(varargin{nextarg}))
      t = varargin{nextarg};
      nextarg = nextarg + 1;
    else
      t = 1:size(sng,2);
    end
  end
  if (length(varargin) >= nextarg)
    options = varargin{nextarg};
  else
    options = struct;
  end
  options = default(options,'sliderwindow',true);
  

  xlim = [min(t) max(t)];
  ylim = [min(f) max(f)];
  [i,j,s] = find(sng);
  s = abs(s);
  snorm = s/min(s);
  lsnorm = log10(snorm);
  spm = sparse(i,j,lsnorm,size(sng,1),size(sng,2));

  if (isfield(options,'WithAmplitude') && options.WithAmplitude)
    hax = SplitVert(0.8);
    hax = hax([2 1]);
    axes(hax(2));
    pow = sum(abs(sng).^2);
    plot(t,pow);
    set(gca,'XLim',xlim,'XTick',[]);
    axupdtfcn = struct('axisupdatefcn',{{@setspimg,@set}});
  else
    hax = gca;
    axupdtfcn = struct('axisupdatefcn',@setspimg);
  end

  cla(hax(1),'reset');
  set(hax(1), 'NextPlot', 'replacechildren', 'TickDir','out','UserData', ...
           struct('m',spm,'xlim',xlim,'ylim',ylim,'climrank',0.9));
  colormap(1-gray);
  axis xy;
  setspimg(hax(1), 'XLim', xlim, 'YLim', ylim);
  
  if options.sliderwindow
    sliderwindow(hax,axupdtfcn);
  end
