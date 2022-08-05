% Analysis of Mouse Videos Using Tracking
%
% Analysis of Mouse Exits/Entrances
%   analyze_tube_traversals.m - extract location information from pixel intensities.
%   track_tube_traversals.m - GUI for quantifying video, entrances and exits to tubes
%   inspect_tube_traversals.m - check movie correspondence to automated behavioral data
%
%
% Analysis of Mouse Center of Mass
%   centerOfMass.m - calc the center of mass of a mask
%   track_mouse_in_ROI_analyze.m - track a single mouse in a ROI
%   track_mouse_in_ROI_gui.m - GUI for behavioral experiments, picking movies and regions
%   track_mouse_in_ROI_inspect.m - check movie correspondence to automated behavioral data
%   track_mouse_in_ROI_twinchambers.m -pick out frames in which mouse is in chamber
%   twinchambers_dist_from_wall.m - compute mouse's distance from each wall in the two-chamber design
%   
%  
%  ? Definitions correct ? Still need editing
%  color2gray.m - change color pixels to grayscale for analysis 
%  detect.m - tracks mouse based on intensity of contrast for each frame 
%  ffmpeg.m - to call mat_ffmpeg in order to jump to frame in video 
%  maskFromInten.m - calc new mask based on maskFromTexture, now use intensity as the hint
%  mask_objs.m - 
%  mmreader_ffmpeg.m - similar to mmreader in matlab
%  movie_quantify_intensity.m - measure RGB intensity in blobs of pixels
%  parse_rois.m - parse roi file saved by roidef
%  test.m -
%  test_thresh.m - to set the contrast for an image
%  to8bit.m - 
%  track_mice.m -  this is a user-friendly wrapper around detect()
% 
%  
%  
%
%
%
%
%
%
%
%
%