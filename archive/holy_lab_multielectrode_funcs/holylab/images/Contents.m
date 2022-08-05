% Functions for image analysis
%
% File IO & format conversion:
%   stackmm             - Fast, convenient access to image stacks
%   imreadheader        - Read image header files
%   imreadheaderOld     - Read old image header files
%   imagine2med         - Convert stacks to various medical file formats
%   stack2movie         - Export stacks as playable movies (AVIs)
%   imdatafilename2headerfilename - Find data file corresponding to
%                                   header
% 
% Display:
%   stackviewer         - A GUI for examining image stacks
%   imsubplot           - Create subplots that don't leave bad gaps
%   imshowrg            - Display two grayscale images in red & green channels
%   colorize_dfof       - Encode deltaF/F with colors
%   imrangegui          - A GUI for choosing color mapping
%
% Filtering and interpolation:
%   imfilter_gaussian   - Quickly filter with a gaussian (multidimensional)
%   imreduce            - Shrink a multidimensional image, after antialiasing
%   imfilter1dnan       - Filter along a single dimension, ignoring NaNs
%   gaussian_filter     - Create a multidimensional gaussian (use
%                           imfilter_gaussian instead)
%   iminterp            - Fast linear interpolation of 1, 2, or 3d images
%
% Processing:
%   image_shift         - Translate an image
%   imdilatend          - Dilate a multidimensional image
%   immean              - Take the mean over a set of images
%   imsvd               - Calculate a space-time decomposition of a movie
%   stack2zcol          - Reshape an image stack
%   zcol2stack          - Reverse the reshaping above
%   
% Image registration:
%   register_multires   - Multiresolution non-rigid registration
%   register_nonrigid   - Optimize registration by pseudo-gradient descent
%   register_warp       - Calculate the registered image
%   register_plotg      - Plot the deformation used to register images
%   register_g0         - Create the identity deformation
%   register_expandg    - Change the scale of a deformation
%     + utility functions: register_jacobian,register_detJ,
%         register_Jtinvw, register_g2detJ. register_graddescent is a
%         slightly older (still highly functional) version of
%         register_nonrigid.
%   register_rigid      - Rigid registration (limited functionality)

%
% High level analysis routines:
%   tissueboundary      - Find the top edge of a piece of tissue
%   imnoiseanalysis     - Analyze properties of cameras & PMTs
%   pyramids_imagine    - Pre-calculate image pyramids for imagine stacks
%
% Simultaneous segmentation & registration (SSR):
%   Needs rehashing with new registration framework
%
% Deconvolution: (of dubious status)
%   deconv_filter_create  - Make a deconvolution filter by Wiener filtering
%   deconv_gaussian     - Approximately invert a gaussian filter
%   deconvstack fcns    - GUIs for deconvolution?
%
% ROI analysis: (someone should write this section sometime)
%
%
% Old File IO & analysis:
%   vimage              - "Virtual image", manipulate images in delayed fashion
%   imload              - Loads real data corresponding to vimages
%   imphys*             - Routines related to the IMPHYS structure
%   imstimcalcdf        - Calculate deltaF/F upon stimulation
%   imstimsvd           - Calculate a space-time decomposition of deltaF
%   imwsum              - Weighted sum of images
%
