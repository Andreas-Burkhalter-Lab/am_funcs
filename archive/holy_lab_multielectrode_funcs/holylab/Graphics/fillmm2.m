function hout = fillmm2(vmin,vmax,varargin)
% FILLMM2: fill region between min and max vectors
% 
% fillmm2(vmin,vmax,x) plots filled polygonal region
% for a function whose envelope is specified by vmin,vmax pairs
% The x coords are given by x.
%
% fillmm2(vmin,vmax) defaults to x = 1:length(vmin)
%
% fillmm2(...,'parameter',value,...) passes on sets of parameter/value
% pairs to the elementary fill function.
%
% h = fillmm2(...) returns the handle to the patch object.
%
% This version works around a bug in the MATLAB fill function.
%
% See also: FILL.

% Have to work around a bug in the MATLAB code, with limit
% on number of points passed to fill
% (oops, this workaround never got finished. Just be sure to decimate
% so that there are fewer than 4000-odd points in x.)

%h = fillmm2work(vmin,vmax,x);
%return;
npts = length(vmin);
if (nargin < 3 | ~isnumeric(varargin{1}))
  x = 1:npts;
  varstart = 1;
else
  x = varargin{1};
  varstart = 2;
end
if (nargin-2 >= varstart)
  pvpairs = varargin(varstart:end);
else
  pvpairs = {};
end
if (length(vmax) ~= npts | length(x) ~= npts)
  error('All inputs must have the same length');
end
imax = 10000;
istart = 1;
counter = 1;
is_hold = ishold;
hold on
while (istart < npts)
  iend = istart+imax;       % There's overlap between last patch & this one
  if (iend > npts)
    iend = npts;
  end
  h(counter) = fillmm2work(vmin(istart:iend),vmax(istart:iend),x(istart:iend),pvpairs);
  counter = counter+1;
  istart = istart+imax;
end
if nargout > 0
  hout = h;
end
if ~is_hold
  hold off
end
return;

function h = fillmm2work(vmin,vmax,x,pvpairs)
yy = [vmin flipdim(vmax,2)];
xx = [x flipdim(x,2)];
if isempty(pvpairs)
  h = fill(xx,yy,'b','EdgeColor','b');
else
  h = fill(xx,yy,'b','EdgeColor','b',pvpairs{:});
end
return
