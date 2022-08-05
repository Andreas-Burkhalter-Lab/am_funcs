function th = legendrgb(hax,labels,colors,options)
% legendrgb: add colored labels to images
% Syntax:
%   texth = legendrgb(hax,labels,colors,options)
% where
%   hax is a handle to the axis
%   labels is a cell array containing the strings you want to show
%   colors is a cell array, one per string
%   options may have the following fields:
%     pstr (default 't'): a "position string" describing where you want the
%       labels. Currently only 't' is supported.
%
% Warning: this function may be a bit incomplete and flaky, and with use it
% seems possible that it will need to change.

% Copyright 2010 by Timothy E. Holy

  if (nargin < 4)
    options = struct;
  end
  options = default(options,'pstr','t',...
    'FontSize',10,'FontWeight','bold');
  th = zeros(1,length(labels));
  xl = get(hax,'XLim');
  yl = get(hax,'YLim');
  n_labels = length(labels);
  switch options.pstr
    case 't'
      if strcmp(get(hax,'YDir'),'normal')
        ybase = yl(2);
      else
        ybase = yl(1);
      end
      for i = 1:n_labels
        th(i) = text(xl(1),ybase,labels{i},'parent',hax,'color',colors{i},...
          'VerticalAlignment','baseline',...
          'FontSize',options.FontSize,...
          'FontWeight',options.FontWeight);
      end
      % Fix the vertical position
      extc = get(th,'Extent');
      if iscell(extc)
        ext = cat(1,extc{:});
      else
        ext = extc;
      end
      top = ext(:,2) + ext(:,4);
      for i = 1:n_labels
        if (top(i) > ybase)
          postmp = ext(i,:);
          postmp(2) = postmp(2)-(top(i)-ybase)/2;
          set(th(i),'Position',postmp(1:2));
        end
      end
      % Space horizontally
      pos = get(th(1),'Position');
      pos(1) = 2*xl(1)-ext(1)+diff(xl)/100;
      set(th(1),'Position',pos);
      if n_labels > 1
        dx = sum(ext(:,3)-ext(:,1));
        gap = (diff(xl)-dx)/(n_labels-1);
        cur = 0;
        for i = 2:n_labels
          cur = cur + ext(i-1,3) + gap;
          postmp = get(th(i),'Position');
          postmp(1) = pos(1)+cur;
          set(th(i),'Position',postmp);
        end
      end
  end
  
      
        