% IMPHYSP: a structure for holding processed imaging physiology data
%
% imphysp structures are collections of variables describing the results
% of processing images.  They collect into a single structure the
% information needed to manipulate images.
%
% The IMPHYSP structure has many of its elements in common with the
% IMPHYS structure, and some of the imphys functions will work on IMPHYSP
% structures.
%
% Relative to the IMPHYS structure, the following are meaningless and
% should be left empty or absent:
%
%   stacktime (timing field)
%   imrange (intensity field)
%
% Relative to the IMPHYS structure, the following are added:
%
% Analysis fields:
%   oimfile: "original image file," the filename containing the original
%     raw (unprocessed) frames/stacks.  Note that xrange, yrange, and
%     zrange refer to coordinates in this file.
%   ostacknum: "original stack numbers," a vector of frame numbers in the
%     original images used to compute the processed frame/stack;
%   ostackweight: the weighting given to each frame/stack in the stacknum
%     vector (e.g., -0.25 might be used when 4 frames are used to compute
%     a background, which is subtracted from 6 "foreground" frames each
%     with weight 0.167);
%   filterwidth: the width (sigma) of the gaussian used to filter the
%     image (empty if no filtering);
%   filtersize: the number of pixels in the gaussian filter;
%   resampmag: the resampling magnification (1 or empty = no resampling);
%   resampmethod: the resampling method used in imresize (default ??);
%
% See also: IMPHYS.
