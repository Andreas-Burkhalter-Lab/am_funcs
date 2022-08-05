function he = barerror(hb,L,U)
% BARERROR: put error bars on a bar graph (even grouped)
% When MATLAB generates a grouped bar graph, inserting error bars
% requires knowing the positioning of the individual bars.  This function
% simplifies that task.
%
% Syntax:
%   he = barerror(hb,E)
%   he = barerror(hb,L,U)
% where
%   hb is the vector of handles for the barseries objects;
%   E is the matrix of errors for the bars (it should have the same size
%     as the Y variable you passed to bar);
%   L and U are lower & upper errorbars;
% and
%   he is the vector of handles to the errorbarseries object.
%
% See also: BAR, ERRORBAR.
  
% Copyright 2006 by Timothy E. Holy
  
  if (nargin < 3)
    U = L;
  end
  for i = 1:length(hb)
    hbp = get(hb(i),'Children');
    xd = get(hbp,'XData');
    X(:,i) = mean(xd([1 3],:))';
    Y(:,i) = get(hb(i),'YData')';
  end
  ish = ishold;
  hold on
  he = errorbar(X,Y,L,U);
  set(he,'LineStyle','none')
  if ~ish
    hold off
  end