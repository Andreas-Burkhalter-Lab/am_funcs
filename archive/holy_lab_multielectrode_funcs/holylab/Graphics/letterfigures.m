function hlabels = letterfigures(hax,labels,offset)
% LETTERFIGURES: attach labels to panels in a figure
% Syntax:
%   hlabels = letterfigures(hax,labels)
%   hlabels = letterfigures(hax,labels,offset)
% where
%   hax is a vector of axis handles;
%   labels{i} is the label used for hax(i);
%   offset is a shift in the default position for label{i}, expressed in
%     normalized figure units; offset(i,:) = [shiftx shifty].

if strcmp(get(hax,'type'),'figure')
  error('Syntax disabled, use (hax,labels) syntax');
end

nlabels = length(labels);
if (nargin < 3 || isempty(offset))
  offset = zeros(nlabels,2);
end
if (size(offset,1) == 1)
  offset = repmat(offset,nlabels,1);
end

% Get axis units for posterity (restore later)
caunits = get(hax,'Units');
set(hax,'Units','normalized');
% Calculate the typical panel size, for obtaining a default shift
axpos = get(hax,'Position');
if length(hax)>1
    axpos = cat(1,axpos{:});
end
defaultoffset = median(axpos(:,3:4),1) .* [-.05,.05];
% Create a hidden axis for writing the labels
% If the parents are multiple figures, we have to create multiple hidden
% axes
hfig = get(hax,'Parent');
if length(hax)>1
    hfig = cell2mat(hfig);
end
[ufig,tmp,parentIndex] = unique(hfig);
for i = 1:length(ufig)
  hlax(i) = axes('Parent',ufig(i),'Units','normalized',...
    'Position',[0 0 1 1],'Visible','off');
end


for i = 1:nlabels
  axpos = get(hax(i),'Position');
  letterpos = [axpos(1),axpos(2)+axpos(4)] + ...
    offset(i,:) + defaultoffset;
  hlabels(i) = text(letterpos(1),letterpos(2),labels{i},...
    'Parent',hlax(parentIndex(i)),...
    'HorizontalAlignment','right',...
    'VerticalAlignment','bottom',...
    'Visible','on',...
    'FontWeight','bold');
end
if length(hax)==1
    set(hfig,'Units',caunits);  % Restore the units setting
else
    set(hax,{'Units'},caunits);
end
