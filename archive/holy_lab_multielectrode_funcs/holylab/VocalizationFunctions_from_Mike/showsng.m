function showsng(varargin)
% showsng: display sonogram
% Calling syntaxes:
%     showsng(filename)
%     showsng(sng,p), where sng and p are loaded from the sng.mat file
%           e.g., "load E8_2sng.mat" will bring these into memory. p may
%           be omitted, but then the axes will not be labelled with
%           physical units, nor is the highlighting available.
%     showsng(...,plotparams) allows one to control how the sonograms are
%       displayed. plotparams is a structure with the following fields:
%          freqband:  controls the range on the frequency axis. If
%            supplied as a scalar, it's interpreted as the maximum
%            frequency displayed. If a 2-vector, then it's interpreted as
%            the band [fmin fmax]. All are supplied in units of kHz.
%          colorscale: a 2-vector giving [min max] values for color
%            saturation, where the log10 is first taken of the input
%            data. If this is not supplied, the full range of the image
%            data is used.
%          fdiff: if true, then the displayed quantity is the
%            finite-difference in frequency of the sonogram rather than
%            the frequency itself (i.e., abs(diff(sng)) rather than sng).
%          minutes: when set true, the time axis is in minutes rather
%            than seconds
%          highlight: when set true, the image will be displayed with
%            "keeper" vocalizations shown in color. These
%            keeper periods are defined by the same criterion used in
%            SONSTATS. If this is set, you must also set the field
%            frange. This mode requires frequency information (so either
%            a filename or the p structure must be supplied).
%          highlighttrange: a 2-vector giving the time range (in seconds)
%            to which highlighting should be restricted. If not
%            specified, the whole time range is used.
%          frange: a 4-vector defining [fcontrol fsignal] (each f is a
%            2-vector of [low high] frequency). This is used if highlight
%            is true.
%          showpower: if true, splits the axis so that a small graph of
%            total power can be shown above. If highlight is set true,
%            a second curve showing the cumulative "keeper" power is also
%            shown in red.
%          labelaxes: if true or absent, causes the axes to be labelled.
%          autoreduce: if true, the image size is reduced to the number
%            of pixels available in the axis (at time of figure
%            creation). This greatly decreases the memory required to
%            represent the image.
%          reducextosize: specifies the desired number of pixels along
%            the x-coordinate. This setting overrides autoreduce in this
%            coordinate, but the y-coordinate will still be auto-reduced.
%
% See also: SONSUM, SONSTATS, LOADMC
  
% Copyright 2001 by Timothy E. Holy <holy@pcg.wustl.edu>

  
filename = [];
if ischar(varargin{1})
  filename = varargin{1};
  load(filename);
elseif isnumeric(varargin{1})
  sng = varargin{1};
end
nextarg = 2;
if (isnumeric(varargin{1}) & nargin > 1 ...
     & isstruct(varargin{2}) & isfield(varargin{2},'scanrate'))
  p = varargin{2};
  nextarg = 3;
end
if (nargin >= nextarg)
  options = varargin{nextarg};
  if ~isfield(options,'minutes')
    options.minutes = 0;
  end
else
  options.minutes = 0;
end
if ~isfield(options,'labelaxes')
  options.labelaxes = 1;
end
% Check now to see if we're finite-differencing, so that
% the frequency coordinates get set at the correct size
if (isfield(options,'fdiff') & options.fdiff)
  %sng = abs(diff(sng));
  %sng = exp(abs(diff(log(sng))));
  lsng = log(sng);
  stepsize = 10;
  sng = exp(abs(lsng(1:end-stepsize,:)-lsng(1+stepsize:end,:)));
end
if exist('p')
  t = linspace(0,p.tacq,size(sng,2));
  f = linspace(0,p.scanrate/2000,p.nfreq);
  unitst = 's';
  if options.minutes
    t = t/60;
    unitst = 'min';
  end
  unitsf = 'kHz';
else
  t = 1:size(sng,2);
  f = 1:size(sng,1);
  unitst = 'arb';
  unitsf = 'arb';
end
newplot;
% Cut off above fmax, if desired
findx = 1:length(f);
if (isfield(options,'freqband') & ~strcmp(unitsf,'arb'))
  if (length(options.freqband) == 1)
    findx = find(f <= options.freqband);
  else
    findx = find(f >= options.freqband(1) & ...
                 f <= options.freqband(2));
  end
end
f = f(findx);
cdata1 = log10(sng(findx,:));
if isfield(options,'colorscale')
  colorscale = options.colorscale;
else
  colorscale = [min(min(cdata1)) max(max(cdata1))];
end
% Instead of using the automatic routines, construct the color
% mapping by hand. That way, the highlighting works more easily.
%imagesc(t,f,lgsng,colorscale);
%axis xy;
%colormap(1-gray)
ihigh = find(cdata1 > colorscale(2));
cdata1(ihigh) = colorscale(2);
ilow = find(cdata1 < colorscale(1));
cdata1(ilow) = colorscale(1);
cdata1 = (colorscale(2)-cdata1)/diff(colorscale);
cdata3 = repmat(cdata1,[1 1 3]);
if isfield(options,'frange')
  bandctl = options.frange(1:2);
  bandsig = options.frange(3:4); 
  fsigi = find(f >= bandsig(1) & f <= bandsig(2));
end 
if (isfield(options,'highlight') & options.highlight)
  % Find regions where the power ratio > 1
  ratsng = SngRatio(sng,p,bandsig,bandctl);
  if isfield(options,'highlighttrange')
    tikill = find(t < options.highlighttrange(1) | ...
                  t > options.highlighttrange(2));
    ratsng(tikill) = 0;
  end
  ti = find(ratsng > 1);
  % Color those pixels red
  cdata3(fsigi,ti,1) = 1;
end
allax = gca;
if (isfield(options,'showpower') & options.showpower)
  mksize = 4;
  hax = SplitVert([0.84 0.85]);
  delete(hax(2)); hax = hax([1 3]);
  axes(hax(1))
  pow = sum(sng(fsigi,:));
  % semilogy(t,pow,'.k','MarkerSize',mksize)
  % plot(t,pow,'.k','MarkerSize',mksize)
  plot(t,cumsum(pow),'k')
  %axis tight
  if (isfield(options,'highlight') & options.highlight)
%    line(t(ti),pow(ti),'LineStyle','none','Marker','.','Color','r', ...
%         'MarkerSize',mksize)
    line(t(ti),cumsum(pow(ti)),'Color','r');
  end
  set(gca,'XTick',[],'XLim',t([1 end]),'Box','off','TickDir','out');
  if (options.labelaxes == 1)
    ylabel('Cum. power');
  end
  slidwincmenu(hax);
  axes(hax(2))
  allax = hax([2 1]);
end
imsize = size(cdata3);
imsize = imsize(1:2);
finalsize = imsize;
if (isfield(options,'autoreduce') & options.autoreduce)
  tmpu = get(gca,'Units');
  set(gca,'Units','pixels');
  pos = get(gca,'Position');
  set(gca,'Units',tmpu);
  finalsize = pos([4 3]);
end
if isfield(options,'reducextosize')
  finalsize(2) = options.reducextosize;
end
if any(finalsize ~= imsize)
  cdata3 = imresize(cdata3,finalsize,'bilinear');
  cdata3 = cdata3(2:end-1,2:end-1,:);  % Cut off parts contaminated by edges
end
himg = image('CData',cdata3,'XData',t([1 end]),'YData',f([1 end]));
slidwincmenu(allax,himg);
set(gca,'TickDir','out','Layer','top')
if (options.labelaxes == 1)
  xlabel(['Time (',unitst,')']);
  ylabel(['Frequency (',unitsf,')']);
end
axis tight
axes(allax(end))    % A hack to make sure titles end up in right place