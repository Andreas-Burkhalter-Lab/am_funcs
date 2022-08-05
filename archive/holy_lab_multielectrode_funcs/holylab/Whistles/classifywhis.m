function [po,lowjumpio,highjumpio,sweepio] = classifywhis(filename)
% CLASSIFYWHIS: calculates whistle parameters and classifies
% Syntax:
%    classifywhis(filename)
%    p = classifywhis(filename)
%    [p,lowjumpi,highjumpi,sweepi] = classifywhis(filename)
% where
%   filename is the name of the .sng file (and it assumes the whistimes
%     file is present;
%   p is the structure whisparams;
%   lowjumpi is the index vector for low jumps;
%   highjumpi is the index vector for high jumps;
%   sweepi is the index vector for simple sweeps.
%
% Use the first form if you want graphical output. The second and third
% form suppress the graphical output.
%
% See also: WHISPARAMS.

% Copyright 2004 by Timothy E. Holy
  [pathname,basename,ext] = fileparts(filename);
  wt = load([pathname,basename,'.whistime'])';
  if isempty(wt)
    if (nargout > 0)
      po = [];
      lowjumpio = [];
      highjumpio = [];
      sweepio = [];
    else
      'Non-whistler'
    end
    return
  end
  [sng,h] = ReadSonogram(filename);
  f = linspace(0,h.scanrate/2,h.nfreq);
  whis = whissnip(sng,h,wt);
%  wp = struct('skipzeros',1,'medfilter',11, 'discardends',1,...
%    'freqrange',[30000 105000],'boxfilt',3);
  wp = struct('skipzeros',1,'medfilter',0, 'discardends',0,...
    'freqrange',[30000 105000]);
  p = whisparams(whis,h,wt,wp);
  %p.filename = filename;
  if (nargout == 1)   % A quick exit
    po = p;
    return;
  end
  % Pick out the ones with significant jumps in frequency
  lowjumpi = find(p.haslj);
  highjumpi = find(p.hashj);
  sweepi = find(~p.haslj & ~p.hashj);
  if (nargout > 0)
    po = p;
    lowjumpio = lowjumpi;
    highjumpio = highjumpi;
    sweepio = sweepi;  
    return;
  end
  % No output is desired, do plotting instead
  %%% commented out the following line because causes error; Matlab doesn't know what 'indxup' is -AM
%   fprintf('%d low jumps (%d up, %d down), %d high jumps, %d simple, %d total\n',length(lowjumpi),length(indxup),length(indxdown),length(highjumpi),length(sweepi),length(p.dt));
  subplot(1,2,1)
  hist(p.maxjump,100);
  %   thresh = input('Choose a threshold: ');
  %   if ~isempty(thresh)
  %     sweepi = find(p.maxjump < thresh);
  %   else
  %     sweepi = 1:length(p.maxjump);
  %   end
  %   lowjumpi = setdiff(p.indx,sweepi);
  %
  psimple = whisparamssubset(p,sweepi);
  subplot(1,2,2)
  hist(p.dtheta(sweepi),100);
  mdt = ceil(max(p.dtheta)/pi);
  set(gca,'XTick',(0:mdt)*pi);
  xtl = cell(mdt+1,1);
  xtl{1} = '0';
  for i = 2:mdt+1
    xtl{i} = [num2str(i-1) '\pi'];
  end
  set(gca,'XTickLabel',xtl);
  set(findobj(gcf,'type','axes'),'TickDir','out');
  set(gcf,'HandleVisibility','off');
  figure
  scatter(cat(2,p.pfs{:}),cat(2,p.dpf{:}),'b.');
  hold on
  scatter(cat(2,p.pfs{:})+cat(2,p.dpf{:}),cat(2,p.dpf{:}),'r.');
  hold off
  figure
  whiscorr(psimple,struct('meanonly',1,'nbins',15,'shuffle',1))
  figure
  for i = 1:length(whis)
    subplot(1,3,1)
    xr = [1 size(whis{i},2)];
    imagesc(xr,f,abs(whis{i}));
    line([xr' xr'],repmat(wp.freqrange,2,1),'Color','r');
    set(gca,'YDir','normal')
    subplot(1,3,2)
    if ~isempty(p.peakfreq{i})
      hline = plot(p.peakfreq{i});
      if ~isempty(find(lowjumpi == i))
        set(hline,'Color','r');
      end
      if ~isempty(find(highjumpi == i))
        set(hline,'Color','g');
      end
    else
      cla
    end
    subplot(1,3,3)
    if ~isempty(p.adjf{i})
      plot(p.adjf{i})
    else
      cla
    end
    fprintf('(%d/%d) %d %g\n',[i length(whis) round(p.maxjump(i)) p.dtheta(i)]);
    pause
  end
  return
  thresh = exp(12);
  sweepi = find(p.vara < thresh);
  lowjumpi = setdiff(p.indx,sweepi);
  lgflag = [0 1 0 0 1 0 0 1];
  whiscorr(p,struct('mode','time','corlen',4,'logflag',lgflag,'meanonly',1))
  psimple = whisparamssubset(p,sweepi);
  figure
  whiscorr(psimple,struct('mode','time','corlen',4,'logflag',lgflag,'meanonly',1))
