% Routines for ray tracing
%
% Demo functions:
%   dospim, spim_pcx, spim_aspheric, check_aspheric
%
% Defining optical components (2-d)
%   pcx:            plano-convex lens (use for cylindrical lenses, too)
%
% Defining and plotting surfaces (2-d)
%   opt2dline:      "planar" surface (d = 2)
%   opt2dcircle:    "spherical" surface
%   opt2daspheric:  aspheric surface
%
% Defining rays:
%   ray
%
% Tracing and analyzing rays:
%   raytrace:       trace rays and plot them on screen
%   raywaist:       measure width of beam waist
%   
% 3-d surfaces (incomplete, largely because of plotting issues):
%   planar
%   spherical
%   cylindrical
%
%
% Utility functions:
%   transmitted_ray: Snell's law for transmitted rays
%   perpproj:       Calculate projection matrix to a plane
%   opt_refrindx:   Calculate refractive index of material at wavelength
%
