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



%   testimqmt           - the test code for multithreaded imqinterp
%   test_noiseanalysis  - A test script for imnoise analysis
%   test_tissueboundary - A script for demonstrating tissueboundary 
%   threshold_stack_gui - THRESHOLD_STACK_GUI M-file for threshold_stack_gui.fig
%   tissueboundary_check - verify that the tissue boundary has been calculated
%   tissueboundary_initialize - initial guess for tissue boundary
%   tissueboundary      - find the top surface of the tissue in fluorescence
%   tissue_surface2d    - find the surface of tissue in 2d images
%   trim_mask_edges     - line a logical array with false "pixels"
%   vimage2struct       - convert the vimage list to pure structures
%   vimagebufferobj     - create new image buffer
%   vimageclear         - clear all vimage data from global memory
%   vimageload          - load data containing vimages
%   vimagesave 			- save data that includes vimages (virtual images)
%   vimbufresize		- resize the vimage buffer
%   vimbufsize			- return the current size of the vimage buffer
%   write_corrected_stacks - utility for balancing intensity, removing stripes, clims, and cropping
%   zcol2stack			- convert an image stack to a z-column matrix
%   zern2mirao1			- the dumb way to optimize actuator voltages
% 						  it will not take into consideration influence functions
%   zernike_values		- compute the Zernike polynomials within a pupil


%   jpm_registration_script  - 
%   jpm_savestack_tif.m
%   killedges			- replaces the edges of an array with a particular value
%   legendrgb 			- add colored labels to images
%   make_coarse_movie   - generate a lower-resolution version of an imagine file
%   make_gaussian_filter - create a spatially-transformed Gaussian
%   make_gaussian_image  - create an image as a sum of Gaussians
%   manipulate_rois     - MANIPULATE_ROIS M-file for manipulate_rois.fig
%						  MANIPULATE_ROIS, by itself, creates a new MANIPULATE_ROIS or raises the existing
%					      singleton*.
%   merge_roi_defs      -
%   movieread			- Virtually read imphys data
%   moviewrite			- write a movie to disk
%   multiple_stack_registration  - this script registers all the stacks
%   optocell_analyze 	- the optical data analog of 'analyze' from ephys land
%   plot_valve_transition - get valve interested
%   prepare_bricks		- divide array (or image) volume into bricks of common size
%   psf_gui				- calculate the PSF by interactively clicking on beads

