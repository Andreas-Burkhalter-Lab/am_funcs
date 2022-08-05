function [indexout,x,y] = concratio(fits,stim1,stim2,options)
% CONCRATIO: compute & plot the effective concentration ratio of two stimuli
%
% Under the assumption of a single ligand present at differing
% concentrations in two different stimuli, compute the concentration
% ratio using the fitting parameters extracted from CONCCOMPARE.
%
% This function basically handles the bookkeeping required to go
% from the technically-optimal but otherwise awkward parameterization of
% the fit in CONCCOMPARE.
%  
% Syntax:
%   [indexout,x,y] = concratio(fits,stim1,stim2,options)
% where
%   fits is the results structure returned by conccompare;
%   stim1 and stim2 are strings containing the name of the two stimuli
%     under consideration (the ratio c1/c2 is computed);
%   options is an optional structure with the following fields:
%     responsive: set to 1 if you want to insist on studying only
%       responsive cells (default: 1);
%     goodfit: set to 1 if you want to insist on studying only
%       cells for which the model is a reasonable fit (default: 1);
%     chifactor: "reasonable fit" is defined as chisq < chifactor*dof
%       (default: 5);
%     strongresp: set to 1 if you want to insist on studying only cells
%       which increase their firing rate by at least some minimal amount
%       (see min_maxr) to at least one stimulus (default: 1);
%     min_maxr: the minimum value allowed for rmax, if strongresp is set
%       to 1 (default: 0.5);
%     nobounds: set to 1 if you want to insist on studying only cells
%       which respond to both stimuli (default: 0);
%     showraw: set to 1 if you want to show the raw
%       concentration ratios (default: 1);
%     showdist: set to 1 if you want to show the
%       distribution of concentration ratios (default: 1);
%     shadedconfintervals: set to 1 if you want to shade the
%       concentration ratio in inverse proportion to its width, to help
%       the narrowest ones stand out best (default: 1);
%     minratio & maxratio: determines the range of the concentration
%       ratio axis on the plot.  More subtly, it acts as a cutoff in
%       cases where one of the two stimuli was ineffective and for which
%       the concentration ratio was given by a bound---see
%       below. Defaults: 0.01 and 100;
%     shadingirange and shadingorange: determines the colorscale used in
%       shading. You'll have to look at the code here, I don't remember
%       what I did at this point.
% and
%   indexout is a vector of indices into the cells of the results
%     structure, indicating which cells passed the various tests and are
%     included in the plots;
%   x is the x coordinate of the distribution of concentration ratios;
%   y is the distribution of concentration ratios.  This is computed by
%     adding up a bunch of gaussians whos widths are given by the error
%     bars on the concentration ratio.  When one of the stimuli is
%     ineffective (and therefore whos concentration ratio is an upper
%     bound), the axis limits are used to help define with width of the
%     gaussian.
%
% See also: CONCCOMPARE, CONCCOMPAREPLOT.

% Copyright 2003 Timothy E. Holy
  
if (nargin < 4)
  options = struct;
end
options = croptions(options);
nexpts = length(fits);
rat = [];
expindx = [];
cellindx = [];
raterrl = [];
raterrr = [];
for i = 1:nexpts
  indx1 = strmatch(stim1,fits(i).uniquelabels,'exact');
  indx2 = strmatch(stim2,fits(i).uniquelabels,'exact');
  if (~isempty(indx1) & ~isempty(indx2)) % if both stimuli tested
    goodindx = 1:size(fits(i).rates,2);
    % If desired, take only cells that respond to at least one stimulus
    if options.responsive
      goodindx = goodindx(find(fits(i).responsive(goodindx)));
    end
    % If desired, take only cells for which the model provides a reasonable fit
    if options.goodfit
      goodindx = goodindx(find(fits(i).chisq(goodindx) < options.chifactor * fits(i).dof));
    end
    % If desired, take only cells with a certain minimal increase in firing rate
    if options.strongresp
      goodindx = goodindx(find(fits(i).rmax(goodindx) > options.min_maxr));
    end
    % Eliminate cases where both are bounds
    bound2 = fits(i).bound([indx1 indx2],goodindx);
    keepers = find(all(bound2,1) == 0);
    goodindx = goodindx(keepers);
    bound2 = bound2(:,keepers);
    % If desired, eliminate cases where only one is a bound
    if options.nobounds
      goodindx = goodindx(find(all(bound2==0,1)));
    end
    if (nargout > 0)
      indexout{i} = goodindx;
    end
    % Compute the concentration ratio
    rat = [rat fits(i).relconc(indx1,goodindx)./fits(i).relconc(indx2,goodindx)];
    % Compute the error bars
    % This is a bit tricky because of the weird parameter arrangements
    for j = 1:length(goodindx)
      % If neither parameter is a bound...
      if all(bound2(:,j) == 0)
        % Pick out the rows of the covariance matrix that apply to this
        % pair. It's possible that one of these was the "best stimulus" and
        % therefore has "no" fitting error.
        boundindx = find(fits(i).bound(:,goodindx(j)));
        notboundindx = 1:1+length(fits(i).uniquelabels);
        for k = length(notboundindx):-1:1
          if ~isempty(strmatch(fits(i).fitlabels{k,goodindx(j)},fits(i).uniquelabels(boundindx)))
            notboundindx(k) = [];
          end
        end
        A = fits(i).invcov(:,:,goodindx(j));
        A = A(notboundindx,notboundindx);
        Clabels = fits(i).fitlabels(notboundindx,goodindx(j));
        C = inv(A);
        Cindx1 = strmatch(stim1,Clabels,'exact');
        Cindx2 = strmatch(stim2,Clabels,'exact');
        err = 0;
        if ~isempty(Cindx1)
          err = C(Cindx1,Cindx1);
        end
        if ~isempty(Cindx2)
          err = err + C(Cindx2,Cindx2);
        end
        if (~isempty(Cindx1) & ~isempty(Cindx2))
          err = err - 2*C(Cindx1,Cindx2);
        end
        raterrl(end+1) = err; % Note this is still on log(rat)
        raterrr(end+1) = err;
      elseif bound2(1,j)
        raterrl(end+1) = Inf;
        raterrr(end+1) = 0;
      else
        raterrl(end+1) = 0;
        raterrr(end+1) = Inf;
      end
      if (raterrl(end) == 0 & raterrr(end) == 0)
        warning('We have a problem')
      end
    end
    % Housekeeping
    expindx = [expindx repmat(i,1,length(goodindx))];
    cellindx = [cellindx goodindx];
  end
end

% Do the analysis for plotting purposes
% (also restricts the domain of "Inf")
ratmin = repmat(options.min_ratio,1,length(rat));
ratmax = repmat(options.max_ratio,1,length(rat));
notinf = find(raterrl < Inf);
ratmin(notinf) = rat(notinf).*exp(-raterrl(notinf));
notinf = find(raterrr < Inf);
ratmax(notinf) = rat(notinf).*exp(raterrr(notinf));
badindx = find(ratmin == 0);
ratmax(badindx) = [];
ratmin(badindx) = [];

newplot
nplots = options.showraw + options.showdist;
cplot = 1;
if options.showraw
  if isfield(options,'axes')
    axes(options.axes(cplot));
  else
    subplot(nplots,1,cplot);
  end
  hline = plot([ratmin; ratmax],[1;1]*(length(ratmin):-1:1));
  if options.shadedconfintervals
    width = log(ratmax) - log(ratmin);
    irange = options.shadingirange;
    windx = find(width < irange(1));
    width(windx) = irange(1);
    windx = find(width > irange(2));
    width(windx) = irange(2);
    colorsat = 1./sqrt(width);
    irange = 1./sqrt(irange);
    slope = diff(options.shadingorange)/diff(irange);
    offset = options.shadingorange(2)-slope*irange(2);
    colorsat = slope*colorsat+offset;
    for kk = 1:length(hline)
      set(hline(kk),'Color',[0 0 0]*colorsat(kk) + [1 1 1]*(1-colorsat(kk)));
    end
    set(gcf,'Color',[1 1 1]);
  end
  set(gca,'XScale','log','XLim',[options.min_ratio options.max_ratio],'YLim',[0 length(ratmin)+1])
  if (nplots > 1)
    set(gca,'XTick',[]);
  end
  set(gca,'YTick',[],'FontSize',14,'Visible','off');
  ylabel('Individual sites','Visible','on')
  cplot = cplot+1;
end
if options.showdist
  if isfield(options,'axes')
    axes(options.axes(cplot));
  else
    subplot(nplots,1,cplot);
  end
  % Add a bunch of gaussians of given widths
  x = logspace(log10(options.min_ratio),log10(options.max_ratio),1000);
  lx = log(x);
  y = zeros(size(x));
  for i = 1:length(ratmin)
    sig = (log(ratmax(i)) - log(ratmin(i)))/2;
    mn = (log(ratmin(i)) + log(ratmin(i)))/2;
    y = y + exp(-(lx-mn).^2/(2*sig^2))/(2*pi*sig);
  end
  semilogx(x,y,'k')
  set(gca,'XLim',[options.min_ratio options.max_ratio]);
  set(gca,'YLim',[0 max(y)]);
  % Shade the most extreme regions
  if options.shadetails
    indx = find(x>10);
    xx = x(indx);
    yy = y(indx);
    xx = [xx xx(end) xx(1)];
    yy = [yy 0 0];
    hpatch = patch(xx,yy,[0.5 0.5 0.5]);
    %set(hpatch,'EdgeColor','none')
    indx = find(x<0.1);
    xx = x(indx);
    yy = y(indx);
    xx = [xx xx(end) xx(1)];
    yy = [yy 0 0];
    hpatch = patch(xx,yy,[0 0 0]);
  end
  set(gca,'FontSize',14,'Box','off')
  xlabel('Effective concentration ratio')
  ylabel('Frequency')
  set(gca,'YTick',[],'TickDir','out','TickLength',[0.02 0.05])
end


%errorbar(1:length(rat),log(rat),raterrl,raterrr)

function options = croptions(options)
if ~isfield(options,'responsive')
  options.responsive = 1;
end
if ~isfield(options,'goodfit')
  options.goodfit = 1;
end
if ~isfield(options,'chifactor')
  options.chifactor = 5;
end
if ~isfield(options,'strongresp')
  options.strongresp = 1;
end
if ~isfield(options,'min_maxr')
  options.min_maxr = 0.5;
end
if ~isfield(options,'nobounds')
  options.nobounds = 0;
end
if ~isfield(options,'showraw')
  options.showraw = 1;
end
if ~isfield(options,'showdist')
  options.showdist = 1;
end
if ~isfield(options,'shadedconfintervals')
  options.shadedconfintervals = 1;
end
if ~isfield(options,'shadingirange')
  options.shadingirange = [0.2 5];
end
if ~isfield(options,'shadingorange')
  options.shadingorange = [1 0.2];
end
if ~isfield(options,'shadetails')
  options.shadetails = 1;
end
if ~isfield(options,'minratio')
  options.minratio = 0.01;
end
if ~isfield(options,'maxratio')
  options.maxratio = 100;
end
