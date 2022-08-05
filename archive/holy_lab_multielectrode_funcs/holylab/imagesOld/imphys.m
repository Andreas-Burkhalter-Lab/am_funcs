% IMPHYS: a structure for holding imaging physiology data
%
% imphys structures are collections of variables describing images.  They
% collect information needed to manipulate images into a single
% structure.  There are several flavors of imphys structures, and some
% functions make use of only subsets of the structure.  The required
% subset should be spelled out in the help for each function.
%
% The fields are described in terms of these subsets.
%
% Dimension fields: (including cropping and transformations)
%   xrange: range ([start end]) of pixels actually present/desired in the
%     image field (can be used to crop on the x-coordinate)
%   x2um: conversion factor between x indices and microns
%   yrange, y2um: similar parameters for y-coordinate
%   zrange, z2um: similar parameters for z-coordinate, if each
%     "image" is in fact a stack.
%   stackz: a vector giving z position of all frames within a stack
%
% Intensity fields:
%   imrange: [minvalue maxvalue] of pixels in the image
%
% Timing fields:
%   stacknum: the frame number (or stack # for 3D)
%   stacktime: the time at which a frame/stack was taken, in seconds,
%     relative to the first stack. (should this be a vector for stacks?
%     or should we introduce a new vector?)
%
% Info fields:
%   tag: a string allowing you to label this item, typically by stimulus
%   date: the date of the experiment
%   stimulus: integer identity of stimulus for this frame/stack.  If
%     stacknum is a vector, then stimulus is also a vector, where each
%     element gives the stimulus identity for each stack
%
% Storage fields (used for retrieving image data):
%   imfile: the file name containing the images
%   headerfile: the file name of auxillary data
%   stackfpos: the file position for the start of the stack (consider
%     making this a vector if we're really doing a stack rather than a frame)
%   imfilefmt: the format of the image file (options: 'raw', 'tif',
%     'stk')
%   width: the width of single images, in pixels
%   height: the height of single images, in pixels
%   depth: the number of frames/stack
%   pixelorder: a vector indicating the order in which pixels are written
%     to disk. For example, if the pixels are written by row, this would
%     be [2 1], whereas pixels written by column would be [1 2]. (correct?)
%   imfilemachfmt: the endian status of the image file (used only for
%     'raw')
%   imfileprec: the precision of the image file (used only for 'raw')
%   camera: the name of the camera 
%   tform: image transform, e.g. for registration (see IMTRANSFORM and
%     CP2TFORM); this is applied before any cropping.
%
% Private fields: (rarely used directly)
%   image: the actual image data (typically, use IMPHYSFETCH to get image data)
%
% As a general rule, "manipulated" images (ones which are not recovered
% by reloading from disk) should not contain the "storage" fields.  See
% IMPHYSCOPY for tools to easily copy subsets of these data.
%
% See also: IMPHYSFROM2D, IMPHYSFETCH, IMTRANSFORM, CP2TFORM.