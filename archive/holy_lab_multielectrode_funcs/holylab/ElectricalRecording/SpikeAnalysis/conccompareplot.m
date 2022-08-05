function conccompareplot(results,options)
% CONCCOMPAREPLOT: display fitting results for different stimuli
% Syntax:
%   conccompareplot(results,options),
% where
%   results is the output structure of CONCCOMPARE;
%   options is an optional structure with the following fields:
%     onepage: if true, all cells shown in a single figure, otherwise hit a
%       key to page through individual cells (default true);
%     showraw: if true, creates a panel with the raw data in one panel
%       (default true);
%     showshift: if true, creates a panel with the data points shifted
%       horizontally to trace out a sigmoid is shown (default true);
%     numbers: a vector of specific cell/channel indices to show (all others
%       ignored); this option supersedes the following options
%     convergeonly: if true, displays only cells for which the minimum was
%       found (default false);
%     responsiveonly: if true, diplays only cells which were judged to be
%       responsive to at least one stimulus (default false);
%
% See also: CONCCOMPARE.

% Copyright 2003 by Timothy E. Holy
if (nargin < 2)
  options = struct;
end
options = ccpoptions(options);

% Pick the units that we will show
if options.convergeonly
    cellindx = find(results.minimumfound);
else
    cellindx = 1:size(results.rates,2);
end
if options.responsiveonly
    cii = find(results.responsive(cellindx));
    cellindx = cellindx(cii);
end
if isfield(options,'numbers')
  cellindx = options.numbers;
end

% Pick the stimuli to show
stimuli = results.uniquelabels;
if ~strcmp(options.stimuli,'all')
  stimuli = options.stimuli;
end

% Layout the figure window
ncells = length(cellindx);
nustim = length(stimuli);
if options.onepage
  nrow = min([8 ncells]);
  ncol = ceil(ncells/nrow);
  nrow = ceil(ncells/ncol);
  nplots = ncells;
else
  nrow = 1;
  ncol = 1;
  nplots = 1;
end
naxes = options.showraw + options.showshift;
ncol = naxes*ncol;

% Do the plotting for each cell
for i = 1:ncells
  showingfit = options.showshift && ~results.fitfailed(cellindx(i));
  if options.showraw
    % Plot the raw data
    if isempty(options.axes)
      if (nrow == 1 && ncol == 2 && ~showingfit)
        ncol = 1;
      end
      subplot(nrow,ncol,(min([i nplots])-1)*naxes+1);
    else
      axes(options.axes((i-1)*naxes+1));
    end
    hold on
    for j = 1:nustim
      indx = strmatch(stimuli{j},results.stimlabel,'exact');
      tmp = errorbarlog(results.stimconc(indx),results.rates(indx,cellindx(i)),results.rateerrs(indx,cellindx(i)));
      set(tmp,'Color',options.colororder(mod(j-1,size(options.colororder,1))+1,:),...
        'LineWidth',options.linewidth(mod(j-1,length(options.linewidth))+1));
      if (isequal(options.linewidth,1) && isfield(results,'beststim') && strcmp(stimuli{j},results.beststim{cellindx(i)}))
        % Highlight the best stimulus with a thicker line
        set(tmp,'LineWidth',2);
      end
      h(j) = tmp(1);
    end
    set(gca,'XScale','log','TickDir','out')
    if (ncells > 1 && isfield(results,'responsive'))
      % Title the axis with the cell # and other information
      titstr = num2str(cellindx(i));
      if results.responsive(cellindx(i))
        titstr = [titstr 'r'];
      end
      if results.minimumfound(cellindx(i))
        titstr = [titstr 'c'];
      end
      title(titstr);
    end
    hold off
    axis tight
  end
  if showingfit
%     if (max(results.relconc(:,cellindx(i))) > 1e6)  % skip obviously-broken ones
%       warning(['Skipping cell ' num2str(i) ', the fit is bad']);
%       continue
%     end
    % Plot the shifted data
    if isempty(options.axes)
      subplot(nrow,ncol,naxes*min([i nplots]));
    else
      axes(options.axes(i*naxes));
    end
    hold on
    % Find factor needed to plot in cOverK units
    bestindx = find(results.relconc(:,cellindx(i)) == 1);
    indx = strmatch(results.uniquelabels{bestindx},results.stimlabel,'exact');
    maxc = max(results.stimconc(indx));
    cOKfac = results.cOverK(cellindx(i))/maxc;
    for j = 1:nustim
      indx = strmatch(stimuli{j},results.stimlabel,'exact');
      indxstim = strmatch(stimuli{j},results.uniquelabels,'exact');
      tmp = errorbarlog(results.stimconc(indx)*results.relconc(indxstim,cellindx(i))*cOKfac,results.rates(indx,cellindx(i)),results.rateerrs(indx,cellindx(i)));
      set(tmp,'Color',options.colororder(mod(j-1,size(options.colororder,1))+1,:),...
        'LineWidth',options.linewidth(mod(j-1,length(options.linewidth))+1));
      delete(tmp(2));
      h(j) = tmp(1);
    end
    axis tight
    if options.showsigmoid
      xlim = get(gca,'XLim');
      x = logspace(log10(xlim(1)),log10(xlim(2)),1000);
      y = results.rmax(cellindx(i)) * (x ./ (1+x));
      if isfield(results,'r0')
        y = y + results.r0(cellindx(i));
      end
      plot(x,y,'k','LineWidth',2);
    end
    set(gca,'XScale','log','TickDir','out')
    hold off
  end
  loc = 'NorthWest';
  if ~showingfit
    loc = 'NorthEastOutside';
  end
  if ~options.onepage
    if options.textoutput
      fprintf('%d  resp: %d (%g)   bound: ',i,results.responsive(cellindx(i)),results.responsivetest(cellindx(i)))
      fprintf('%d ',results.bound(:,cellindx(i)))
      fprintf('fit: %g (%d dof)\n',results.chisq(cellindx(i)),results.dof)
    end
    if options.showkey
      legend(h,stimuli{:},'Location',loc);
    end
    if (i < ncells)
      pause
      clf
    end
  end
end
if (options.onepage & options.showkey)
  legend(h,stimuli{:},'Location',loc);
end


function options = ccpoptions(options)
if ~isfield(options,'onepage')
  options.onepage = 1;
end
if ~isfield(options,'showraw')
    options.showraw = 1;
end
if ~isfield(options,'showshift')
  options.showshift = 1;
end
if ~isfield(options,'showkey')
    options.showkey = 1;
end
if ~isfield(options,'showsigmoid')
  options.showsigmoid = 0;
end
if ~isfield(options,'convergeonly')
    options.convergeonly = 0;
end
if ~isfield(options,'responsiveonly')
    options.responsiveonly = 0;
end
if ~isfield(options,'stimuli')
  options.stimuli = 'all';
end
if ~isfield(options,'axes')
  options.axes = [];
end
if ~isfield(options,'textoutput')
  options.textoutput = 0;
end
if ~isfield(options,'colororder')
  options.colororder = get(gca,'ColorOrder');
end
if ~isfield(options,'linewidth')
  options.linewidth = 1;
end