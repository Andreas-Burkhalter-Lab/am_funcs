% MSAMS_VAL_GRADIENT: compute penalty function for MSAMS neighborhoods
%
% MSAMS provides a statistical criterion for the definition of the
% "neighborhood" of a point.  It does so in terms of a penalty function,
% the mean-to-s.e.m. ratio of the size of a mean shift step.  When the
% mean-to-s.e.m. ratio significantly exceeds 1, that is a
% statistically-significant indication that the data show a local trend.
%
% This function does not directly return the mean-to-s.e.m. ratio,
% because that can easily be computed from two other quantities, the
% "step size" and "variance", which this function does return.  The step
% size is also interesting in its own right.
%
% Syntax:
%      [step,stepvar,weights,gradE] = msams_val_gradient(x,x0,kvec,max_threads)
% where
%   x is a d-by-N matrix of data points (one in each column);
%   x0 is the "probe point" (or seed location), a d-by-1 column vector;
%   kvec contains the coefficients of the metric, where the square
%     distance d2 between any two points x and y is expressed as
%           d2(x,y) = sum(kvec.^2 .* (x - y).^2)
% and
%   step contains the components of the "step size" that would be taken
%     from x0, given kvec (a d-by-1 vector)
%   stepvar contains the coordinate components of the uncertainty in step
%     size (a d-by-1 vector)
%   w contains the weight assigned to each point,
%           w_i = exp(-d2(x,y))
%   gradE is the gradient of the penalty function, which can be written
%                 sum(kvec.^2 .* step.^2)
%          E =   -------------------------
%                 sum(kvec.^2 .* stepvar)
%
% From the output, you can thus readily compute E(kvec|x0).  A
% minimally-significant step would be one for which E transitions from
% being < 1 to > 1 as |kvec| is decreased.  (But you want to make sure
% that this transition isn't just noise, so you probably want to insist
% that it rises to something like 10 and then backtrace to 1; however, E
% tends to be much smoother than the sharp-cutoff version, and thus this
% precaution may not be necessary.)
%
% For speed, this is implemented as a multithreaded MEX file.

% Copyright 2008 by Timothy E. Holy
