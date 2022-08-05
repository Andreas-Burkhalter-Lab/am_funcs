% CAMERADATA: a structure for specifying camera performance
% The value I (intensity) reported by the camera for a particular pixel is
% related to the physical measurement in the following way:
%     I = gain * n + eta + bias,
% where n is the # of electrons in the well, gain is the gain factor (digital
% #s/electron), eta is a noise term (readout + quantization noise), and
% bias is the offset bias (the average value reported in the dark). This
% model ignores "bleed-through," patterned noise, and other types of
% artifacts; it basically corresponds to an ideal camera model.
%   
% This structure has the following fields:
%   gain
%   bias
%   sigmap: the standard deviation of eta. Notice this is measured in units
%     of digital #s, not in units of electrons! Under real circumstances it
%     will probably always be >= 1, because of quantization noise.
%   goodpixels (optional): a logical array the size of the image produced by the
%     camera, with value true for pixels that are working and false for
%     dead pixels.
%   pixel_spacing (optional): the physical dimensions of a pixel, in units
%     of your sample dimensions (i.e., including the effect of the optical
%     system). If you're taking an image stack, the 3rd dimension of this
%     can be used to encode the z-spacing between images in the stack.
%   QE (optional): the quantum efficiency of the camera (not needed in most
%     cases).
%
% This structure can be accepted as an input by a number of routines.
%
% See also: IMWIENERDECONV.
