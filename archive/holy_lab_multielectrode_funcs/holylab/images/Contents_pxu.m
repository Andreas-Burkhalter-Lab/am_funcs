% Functions for image analysis
%
% File IO & format conversion:
%   STACKMM			-memory-mapped image stack access class
%   stackheaderfilename2camfilename
%   STACK2MOVIE			-create AVI movies from a stackmm object

% 
% Display:
%   STACKVIEWER			-viewing tool for 4D image stacks
%   stackviewer_wm		-volumetric image data browser using the wheel mouse
%   STACK_STIM_RGB		-overlay stimulus responses in an rgb image
%
% Filtering and interpolation:

%
% Processing:
%   STACK2ZCOL			-convert an image stack to a z-column matrix
%   SORT_CELLS_INTO_REGIONS 	-tests whether each cell_point lies within the polygonal regions
%   rotatestack.m		-rotate a 3-d stack of images
%   ROTATESTACK_CROPGUI		-a tool for assisting choice of crop rectangles
%   ROTATESTACK_PREPARE_IMAGINE -set up data for stack rotation with IMAGINE files
%   
% Image registration:
%   register_frames             - 2D registration within a stack
%   STACK_RESGISTER_TRANSLATION -translation-registers a stackmm object with given options
%   SSR				-simultaneous segmentation and registration

% Stimulus
%   DFOFBYSVD			-compute scalar response using SVD
%   FIND_STIMULUS_START     	-locate the first frame, organized by stimulus presentation
%   GET_TRIALS_DFOF		-a simple function to get temporal traces from imaging data
%   GET_TRIALS_DFOF_CORRECTED   -snip temporal traces, subtracting a control timecourse
%   intensity_corrected         -subtract control responses from intensities vs. time
%
% High level analysis routines:

%
% Simultaneous segmentation & registration (SSR):

%
% Deconvolution: (of dubious status)

%
% ROI analysis: (someone should write this section sometime)
%   reorient_roidef_file	-add a new field into roi_def file for true (tissue) coordinates
%   roi_by_imflow.m		-creates regions of interest in a .imagine movie based on activity
%   roi_crawl			-GUI to view ROI plots, uses roi_traceviewer script
%   roi_merge 			-allows interactive merging of ROIs selected by roi_by_imflow
%   roi_traceviewer             -plots of the roi intensity matrix calculated by calc_stack_roi_inten
%   roi_traceviewer_2
%   roiapply			-measure ROI intensities across entire files
%   roidraw			-draw regions of interest in an axis
%   roimeasure			-measure ROI intensity in an image or sequence of images
%   roiplot			-display ROIs
%   roisturct			-a structure for holding ROI information
%   set_roi
%   single_stack_roi		-manually place ROIs in a single stack
%   SPLIT_ROI_DEFS		-convert ROI structure to cell array of single ROIs
%   stack_roi                   -gui to click roi, grey version
%   stack_roi_rg                -gui to click roi, grey/rg version
%   stack_vol_roi_rg		-gui to click roi, 3D version






  
                                