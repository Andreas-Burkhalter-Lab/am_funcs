% IMPHYSPLOTPARAMS: A structure to specify imphys plot parameters
%
% You specify the type of plot through the field 'mode', which can
% take the following values:
%   'image':         show the (enhanced) raw image
%   'image w/ roi':  show the (enhanced) raw image, with ROIs overlaid
%   'stimulus':      show the stimulus as a function of time
%   'psth':          show the response amplitude in ROIs as a function of
%                    time
%
% Each of these different modes takes auxillary parameters specified as
% additional fields, which will be discussed below.
%
% Image mode:
%   bgindx (optional): an index into the imphys structure array giving
%     the frames used in computing the background (if not present,
%     background subtraction is not performed);
%   fgindx: an index into the imphys structure array giving the frames
%     used in computing the "foreground," i.e. the positive portion of
%     the response;
%   cropstruct (optional): region of the image to show (default all), see
%     IMPHYSCROP for the structure;
%   spatialfilter (optional): filter the image using IMFILTER; for a gaussian
%     filter, you may supply a scalar 'sigma' (see IMPHYSMANIP
%     for more detail) or more generally a 2-d array (see FSPECIAL for
%     some good candidates);
%   spatialresample (optional): set a magnification
%   clim (optional): the intensity limits used in plotting with IMAGESC.
%
% In image mode, the return parameter ("ret") of IMPHYSPLOT is the image
% data.
%
%
% Image w/ ROI mode:
% The parameters are the union of the 'image' mode parameters and the
% 'psth' mode parameters.  The return parameter ("ret") is the image
% data.
% 
%
% Stimulus mode:
%   type (default 'bar'): plotting style used to show stimulus
%     timing. 'bar' produces colored bars which span the duration of the
%     valve opening.
%   toffset (default 0): value added to all frame times, can be used to
%     adjust the time origin.
%   valvecolor (optional): a structure array, with fields 'valvenum'
%     (giving the valve number) and 'color' (giving the corresponding
%     color, as a string or as an RGB triple).
%
% In stimulus mode, the return parameter ("ret") of IMPHYSPLOT is a
% vector of handles to the plot objects.
%
%
% PSTH mode:
%   roidata: a structure array specifying the geometrical data for each
%     ROI (see PLOTROI);
%   roioptions: an optional structure array with information about how
%     ROIs should be plotted (see PLOTROI).
