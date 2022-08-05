function [hh]=spsngplot(sng,varargin)
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
%     trange: specify as [start stop] (in seconds) for the time range to
%       show. If supplied, sliderwindow defaults to false.
%     hfig: if supplied, plots in the chosen figure (supply as handle to
%       figure)
%     hax: if supplied, plots in the chosen axis (supply as a handle to the
%       axis)
%     scale:  absolute scale you specify instead of a normalized one

% Copyright 2002 by Timothy E. Holy <holy@pcg.wustl.edu>

  
  % Argument parsing and defaults
  if ischar(sng)
    filename = sng;
    [sng,header] = ReadSonogram(filename);
    f = linspace(0,header.scanrate/2000,header.nfreq);
    t = linspace(0,header.nscans/header.scanrate,header.columnTotal);
    nextarg = 1;
  else
    nextarg = 1;
    if (length(varargin) >= nextarg && ~isstruct(varargin{nextarg}))
      f = varargin{nextarg};
      nextarg = nextarg + 1;
    else
      f = 1:size(sng,1);
    end
    if (length(varargin) >= nextarg && ~isstruct(varargin{nextarg}))
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
  options = default(options,'sliderwindow',~isfield(options,'trange'));

  xlim = [min(t) max(t)];
  ylim = [min(f) max(f)];
  [i,j,s] = find(sng);
  s = abs(s);
  if ~isfield(options,'scale')
    snorm = s/min(s);
    lsnorm = log10(snorm);
  else
    s(s>options.scale(2))=options.scale(2);
    s(s<options.scale(1))=options.scale(1);
    lsnorm=s;
  end

  spm = sparse(i,j,lsnorm,size(sng,1),size(sng,2));

  if isfield(options,'trange')
    xlim=options.trange;
    pa = find(t<options.trange(2));
    rt = find(t>options.trange(1));
    part = intersect(pa,rt);
    spm=spm(:,part);
  end
  
  if isfield(options,'hax')
    hax = options.hax;
  else
    if isfield(options,'hfig')
      figure(options.hfig);
    end
    hax = gca;
  end
  
  if (isfield(options,'WithAmplitude') && options.WithAmplitude)
    hax = SplitVert(0.8,hax);
    hax = hax([2 1]);
    pow = sum(abs(sng).^2);
    plot(hax(2),t,pow);
    set(hax,'XLim',xlim,'XTick',[]);
    axupdtfcn = struct('axisupdatefcn',{{@setspimg,@set}});
  else
    axupdtfcn = struct('axisupdatefcn',@setspimg);
  end

  cla(hax(1),'reset');
  if ~isfield(options,'scale')
    set(hax(1), 'NextPlot', 'replacechildren', 'TickDir','out','UserData', ...
      struct('m',spm,'xlim',xlim,'ylim',ylim,'climrank',0.9));
  else
    set(hax(1), 'NextPlot', 'replacechildren', 'TickDir','out','UserData', ...
      struct('m',spm,'xlim',xlim,'ylim',ylim,'clim',options.scale));
  end
  
  colormap(1-gray);
  axis xy;
  setspimg(hax(1), 'XLim', xlim, 'YLim', ylim);
  
  if options.sliderwindow
    sliderwindow(hax,axupdtfcn);
  end
