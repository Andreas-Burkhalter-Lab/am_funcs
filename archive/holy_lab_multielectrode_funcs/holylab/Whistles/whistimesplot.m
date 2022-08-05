function whistimesplot(sng,header,twhis,featurevecs)
% WHISTIMESPLOT: plot whistle times and features
% Syntax:
%   whistimesplot(sng,header,twhis,featurevecs)
% where
%   sng, header come from READSONOGRAM
%   twhis is a 2-by-n vector giving start and end whistle times (see
%     WHISTIMES);
%   featurevecs (optional) is a structure array with fields 'name',
%     'value', 'thresh' containing the features computed from the
%     sonogram, usually of interest for identifying whistles.
%
% See also: WHISTIMES.
  
% Copyright 2004 Timothy E. Holy.
  
  f = linspace(0,header.scanrate/2,header.nfreq);
  t = linspace(0,header.nscans/header.scanrate,header.columnTotal);

  nplots = 2;
  features = 0;
  if (nargin > 3)
    features = 1;
    nplots = nplots+length(featurevecs);
  end
  hax(1) = subplot(nplots,1,1);
  xlim = [min(t) max(t)];
  ylim = [min(f) max(f)];
  [i,j,s] = find(sng);
  spm = sparse(i,j,log(abs(s)/header.threshold),size(sng,1),size(sng,2));
  set(hax(1), 'NextPlot', 'replacechildren', 'UserData', ...
              struct('m',spm,'xlim',xlim,'ylim',ylim,'climrank',0.9));
  colormap(1-gray);
  axis xy;
  setspimg(hax(1), 'XLim', xlim, 'YLim', ylim);
  if (features)
    for i = 1:length(featurevecs)
      hax(i+1) = subplot(nplots,1,i+1);
      plot(t,featurevecs(i).value);
      axis tight
      ylabel(featurevecs(i).name)
      line(xlim,featurevecs(i).thresh * [1 1],'Color',0.4 * [1 1 1]);
    end
  end
  hax(nplots) = subplot(nplots,1,nplots);
  cla
  sz = size(twhis);
  if (sz(2) == 2 & sz(1) > 2)
    twhis = twhis';
  end
  plot(twhis, ones(size(twhis)), 'r');
  set(gca,'XLim',xlim)
  updtfnc{1} = [@setspimg];
  for i = 2:nplots
    updtfnc{i} = @set;
  end
  sliderwindow(hax,struct('axisupdatefcn',{updtfnc}));
