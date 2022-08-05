% Functions for image analysis
%
% File IO & format conversion:
%   stackmm             - Fast, convenient access to image stacks
%   IMAGINE2MED         - Convert imagine files to medical image file formats (IMG, DICOM, NII, DCM)
%   Imagine2multitiff   - Take .imagine binary data and save it down in multi-tiff format
%   imdatafilename2headerfilename - Find data file corresponding to
%                                   header
%   IMFILE              - prepare to read image data from a file
%   imreadheader        - Read image header files
%   imreadheaderOld     - Read old image header files
%   imagine2med         - Convert stacks to various medical file formats
%   Import_zstack_from_numbered_tifs   - Allows importing into image structure for future analysis
%   stack2movie         - Export stacks as playable movies (AVIs)
%   Imsaveas            - Work around MATLAB bug in SAVEAS
%   IMWRITE             - Save an image in raw format
%   imdatafilename2headerfilename - Find data file corresponding to header
                               
% Display:
%   stackviewer         - A GUI for examining image stacks
%   imsubplot           - Create subplots that don't leave bad gaps
%   imshowrg            - Display two grayscale images in red & green channels
%   IMSHOWRGB           - Display an overlay of images in red/green/blue channels
%   colorize_dfof       - Encode deltaF/F with colors
%   imrangegui          - A GUI for choosing color mapping 
%   IMSHOWSC            - Show scaled grayscale image  
%   IMSUBPLOT           - Create axes for image subplots using screen space efficiently

% Filtering and interpolation:
%   imfilter_gaussian   - Quickly filter with a gaussian (multidimensional)
%   imreduce            - Shrink a multidimensional image, after antialiasing
%   imfilter1dnan       - Filter along a single dimension, ignoring NaNs
%   gaussian_filter     - Create a multidimensional gaussian (use
%                           imfilter_gaussian instead)
%   iminterp            - Fast linear interpolation of 1, 2, or 3d images
%   IMQINTERP           - Quadratic interpolation and gradients for images

% Processing:
%   image_expand        - Takes an image and expands the exact image to a larger size
%   image_shift         - Translate an image
%   image_snip_qinterp  - Cut out rectangular regions of an image, with subpixel interpolation
%   IMCOMPARE           - Consistently evaluate differences between images(calculates the 
%                         mean-square-error between a reference image and one or more test images)
%   IMDIALEND           - Multidimensional single-pixel dilation of images 
%   imdilatend          - Dilate a multidimensional image
%   immean              - Take the mean over a set of images
%   imsvd               - Calculate a space-time decomposition of a movie
%   stack2zcol          - Reshape an image stack
%   zcol2stack          - Reverse the reshaping above
%   IMRS                - Image resample at lower resolution
%   IMMEANSHIFT         - Shift points to the local center of mass of an image
%   IMWSUM              - Calculate the weighted sum of images
%   IMSTIMCALCDF        - Image fluorescence changes (deltaf) upon stimulation
%   IMSTIMCALCSVD       - Analyze stimulus responses by SVD
%   IMSTIMVIEWDF        - Show and analyze deltaf changes upon stimulation
%   IMSVD               - Analyze images by SVD
%   IMWIENERDECONV      - Wiener deconvolution of images with noise handling
%   IM_GAUSSIANS_ERR    - Compute the fitting error between objects and Gaussians,
%                         error measure includes all sufficiently-nonzero model pixels
%   im2gaussians        - Find gaussian-shaped objects in an image

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
%   IMNOISEANALYSIS     - Analyze properties of cameras and PMTs (gain, bias and noise)

% Simultaneous segmentation & registration (SSR):
%   Needs rehashing with new registration framework
%
% Deconvolution: (of dubious status)
%   deconv_filter_create  - Make a deconvolution filter by Wiener filtering
%   deconv_gaussian     - Approximately invert a gaussian filter
%   deconvstack fcns    - GUIs for deconvolution?
%
% ROI analysis: (someone should write this section sometime)
%   intensity_stepper   - Opens a keystepper-controlled window of your ROIs

%
%
% Old File IO & analysis:
%   imload              - Loads real data corresponding to vimages
%   imphys*             - Routines related to the IMPHYS structure
%   imstimcalcdf        - Calculate deltaF/F upon stimulation
%   imstimsvd           - Calculate a space-time decomposition of deltaF
%   imwsum              - Weighted sum of images


%   ????
%   imagine_stim_lookup - Finds valve numbers and labels for a .imagine file
%   imshowrg_stim_gui   - In case imagine crashed in the middle of acquisition, you 
%                         have to go and  change the .imagine file's stack number field 
%                         to represent the last working stack.

%   imqinterp_setthread_num
%   interpolate3dTform
%   intensity_batch

