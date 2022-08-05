function varargout = addbreaks(varargin)
% ADDBREAKS: help plot lines with breaks at missing values
% The x axis points must have even spacing, except for the possibility of
% a few missing values.
% Syntax:
%   [xo,yo] = addbreaks(x,y);
%   [xo,yo,lo,uo] = addbreaks(x,y,l,u);
%   [xo,yo,...] = addbreaks('break',breakindx,x,y,...);
% where
%   x is a vector of positions;
%   y,l,u,... are vectors or matrices (columnwise) with the same numbers
%     of points;
% and
%   xo,yo,... are vectors/matrices with NaNs inserted at points where
%     the x-coordinate changes by anything other than the minimum
%     inter-point spacing.  Note that breaks will automatically be
%     introduced at places where the x-coordinate decreases.
%
% Alternatively, breaks can be specified to occur before particular
% points, labelled by breakindx.
%
% The output points can be used in LINE, PLOT, or ERRORBAR to create
% broken-line plots.
  
% Copyright 2005 by Timothy E. Holy
  
  argoffset = 0;
  % Find breaks
  if (ischar(varargin{1}) && strcmp(lower(varargin{1}),'break'))
    ibreak = [1; varargin{2}(:); Inf];
    argoffset = 2;
  else
    x = varargin{1}(:);
    dx = diff(x);
    dx = min(dx(dx > 0));           % To allow backward jumps in value, too
    xint = round((x-x(1))/dx);
    dxint = diff(xint);
    ibreak = [1; find(dxint > 1)+1; length(x)+1];
  end
  
  nlines = length(ibreak)-1;
  % Insert NaNs
  for i = 1:length(varargin)-argoffset
    tmp = varargin{argoffset+i};
    sztmp = size(tmp);                % Save original shape info
    if isvector(tmp)
      tmp = tmp(:);
    end
    varargout{i} = nan(size(tmp,1)+nlines-1,size(tmp,2));
    for j = 1:nlines
      varargout{i}(ibreak(j)+j-1:ibreak(j+1)+j-2,:) = ...
          tmp(ibreak(j):ibreak(j+1)-1,:);  % Copy data
    end
    if (isvector(tmp) && sztmp(2)>sztmp(1))
      varargout{i} = varargout{i}';   % Restore original shape
    end
  end