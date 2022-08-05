function sonsum(filenames,titles,plotparams)
% sonsum: summarize sonograms for several files
% totpow = sonsum(filenames,titles,plotparams) 
%   where filenames (optional) is a cell array of *sng.mat files (if
%     empty, uses UIGetFiles GUI)
%   titles (optional) is a cell array of string titles for each graph
%   plotparams (optional) controls the way the plot is drawn. See SHOWSNG.
%
% See also: SHOWSNG.
  
% Copyright 2001 by Timothy E. Holy <holy@pcg.wustl.edu>

if (nargin == 0 | isempty(filenames))
  filenames = UIGetFiles;
end
if ~iscell(filenames)
  filenames = {filenames};
end
nfiles = length(filenames);
if (nargin < 2 | isempty(titles))
  titles = filenames;
end
% Supply default plotting parameters
if (nargin < 3 | ~isfield(plotparams,'colorscale'))
  % multi-axis plots should always be on a consistent colorscale
  plotparams.colorscale = [-4 -2];  
end
if ~isfield(plotparams,'frange')
  plotparams.frange = [25 35 55 110];
end
if ~isfield(plotparams,'highlight')
  plotparams.highlight = 1;
end
if ~isfield(plotparams,'showpower')
  plotparams.showpower = 1;
end
plotparams.labelaxes = 0;
bandsig = plotparams.frange(3:4);
bandctl = plotparams.frange(1:2);
clf
% Set up axis labels
% Use trick of invisible axes so later
% subplot commands won't delete the axis labels
hax1 = gca;
hy1 = ylabel('Frequency (kHz)');
hx1 = xlabel('Time (s)');
set(hax1,'Visible','off','HandleVisibility','off');
set([hx1 hy1],'Visible','on','HandleVisibility','on');
% Now get down to business
% Check to see if we can skip x-axis ticks
% (if all files are of the same duration)
for i = 1:nfiles
  load(filenames{i},'p');
  tacq(i) = p.tacq;
end
xticks = 1;
if all(tacq == tacq(1))
  xticks = 0;
end
for i = 1:nfiles
  subplot(nfiles,1,i);
  showsng(filenames{i},plotparams);
  htit = title(titles{i});
  set(htit,'Interpreter','none');
  set(gca,'Tag','son','TickDir','out');
  if (i < nfiles & ~xticks)
    set(gca,'XTick',[]);
  end
end
set(hax1,'HandleVisibility','on');
