% REGISTER_WARP_UTILITIES: do warping computations of image registration
%
% This is a MEX function that can handle a variety of computations that
% occur in image registration.  Eventually it will be multithreaded for
% enhanced performance on SMP machines.  These multiple computations are
% arranged in a single function because there is a fair amount of overlap in
% the tasks, and this seems to be a good way of both preventing code
% divergence and reducing redundant computation.
%
% General syntax:
%   [output1,...] = register_warp_utilities(data)
% where data is a structure.  The data types, and their values, control
% the calculation performed by this function.
%
% Here are the types of calls you can make:
% Calculating the warped image (obtained with data.output = 'w'):
%   img = register_warp_utilities(data)
% where
%   data.psim is the "moving" image (or sqrt(moving image), see below),
%     must be of type single;
%   data.g is the deformation (see REGISTER_G0 for syntax);
%   data.covariant is boolean, if true it uses covariant deformation;
%   data.sqrt is boolean, if false it treats psim (and psit if present,
%     see below) as a real image, if true it treats psim as
%     sqrt(real image).
%
% Calculating the fitting error and the gradient of the fitting error
% (obtained with data.output = 'eg'):
%   [err,graderr] = register_warp_utilities(data)
% where all the fields listed above are possible, as well as
%   data.psit is the "target" image (or sqrt(target image) as above);
%   data.gradpsit (optional) is the (spatial) gradient of psit (you can
%     save some recomputation by supplying this);
%   data.sigma (optional) if present, will perform smoothing on the
%     gradient before returning (sigma should be a vector, with one entry
%     per dimension, giving the amount of smoothing; see
%     IMFILTER_GAUSSIAN).
% With this syntax, you can also get just the error by asking for only
% the first output, saving some time.
%
% Setting data.output = 'egw' gives the warped image as the third output,
% and so on.
%
% See also: REGISTER_G0, IMFILTER_GAUSSIAN.

% Copyright 2006 by Timothy E. Holy
