==> Contents.m <==
% Routines for ray tracing

==> DMmanipulator.m <==
% DMmanipulator: takes old points of the deformable mirror and a point
% where the rays have to converge, then deforms the mirror to fit that

==> aspher_eq.m <==
% standard equation for an aspheric lens surface

==> calculate_aberrations.m <==
% calculate the classical aberrations of a set of rays

==> check_aspheric.m <==
% No introduction so far

==> create_inffoc_rays_for_aberrations.m <==
% create rays for testing infinity-focused optics

==> cylindrical.m <==
% No introduction so far

==> defoc_coefs_planar.m <==
% defocus coefficients for planar "tissue"

==> doao.m <==
% script to study affect of single deformable mirror in adaptive optics

==> dospim.m <==
% No introduction so far

==> generate_rays.m <==
% creates a collection of rays that passes through a lens

==> grindebug.m <==
% No introduction so far

==> grin_calculations.pdf <==
% cannot be opened

==> lightsheet_from_collimated.m <==
% A script for testing different lightsheet configurations, assuming we
% have as input a collimated beam (hopefully achromatized)

==> lightsheet_opt.m <==
% No introduction so far

==> lightsheet_pcx.m <==
% light sheet illumination using plano-convex collimator

==> lightsheet_traceconfig.m <==
% do raytracing through lightsheet illumination

==> linearlens.m <==
% an ideal (sine-condition satisfying) lens model
% This implements the linear model in the Holy manuscript.  For the
% linear model, this is more flexible than RAYTRACE_IDEAL_LENS, which
% assumes the optic axis is along the z-axis.

==> objective_aspheric.m <==
% No introduction so far

==> objective_aspheric_trace.m <==
% No introduction so far

==> objective_grin.m <==
% No introduction so far

==> objective_grin_trace.m <==
% No introduction so far

==> opt2dDM.m <==
% OPT2DMIRROR: defines a deformable mirror for ray tracing

==> opt2dDMflat.m <==
% OPT2DDMFLAT: a simple deformable mirror for ray tracing

==> opt2dDMspline.m <==
% No introduction so far

==> opt2daspheric.m <==
% OPT2DASPHERIC: define an aspheric "surface" (line) for ray tracing

==> opt2dcircle.m <==
% OPT2DLINE: define a circular "surface" (line) for ray tracing

==> opt2dgrin.m <==
% OPT2DGRIN: define a grin lens for ray tracing

==> opt2dline.m <==
% OPT2DLINE: define a planar "surface" (line) for ray tracing

==> opt2dmirror.m <==
% OPT2DMIRROR: define a flat mirror for ray tracing

==> opt2dparabolicmirror.m <==
% No introduction so far

==> opt2dperfect.m <==
% paraxial approximation for all rays (even non-paraxial)

==> opt2dprojective.m <==
% paraxial lens approximation for all rays

==> opt2dsinecond.m <==
% trace rays as required by the sine condition

==> opt2dspline.m <==
% define an arbitrary interface for ray tracing

==> opt3d_raytrace_spherical_onaxis.m <==
% 3d-raytracing for on-axis spherical surfaces

==> opt_refrindx.m <==
% refractive index calculations.

==> opticalAperture.m <==
%  location, orientation, and aperture information for optical components
% This is the base class for all optical components. It handles the
% geometry and transformations as well as tracing of rays to the aperture
% plane. 

==> opticalQuadratic.m <==
% opticalQuadratic: spherical, cylindrical, and planar optical surfaces

==> opticalRayBundle.m <==
% a collection of rays in a common coordinate system.
% This defines rays via a position, direction of propagation, wavelength,
% and intensity.

==> optical_collimated_bundle.m <==
% build a collimated beam of rays

==> optical_flip_component.m <==
% flip a component along its optic axis

==> optical_layout_onaxis.m <==
% building lens systems with a single optic axis

==> optical_plot_components_2d.m <==
% draw two-dimensional projections of components

==> optical_plot_components_3d.m <==
% draw component surfaces in three dimensions

==> optical_plot_rays_2d.m <==
% plot traced rays in a plane containing the optic axis

==> optical_plot_rays_3d.m <==
% optical_plot_rays_3d: plot traced rays

==> optical_point_bundle.m <==
% optical_point_bundle: build a beam of rays emerging from a single point

==> optical_rotate_component.m <==
% optical_rotate_component: rotate the optic axis of a component

==> optical_shift_component.m <==
% optical_shift_component: change the position of a component

==> optical_stock_lenses.m <==
% optical_stock_lenses: find and use stock lenses

==> optical_trace.m <==
% optical_trace: ray tracing through component layouts

==> optical_trace_onaxis.m <==
% optical_trace_onaxis: ray tracing for layouts with a single optic axis

==> parametrize_inffoc_lens.m <==
% No introduction so far

==> pcx.m <==
% PCX: make a plano-convex lens

==> perpproj.m <==
function P = perpproj(v)
% PERPPROJ: calculate the matrix to project perpendicular to a vector
% Syntax:
%   P = perpproj(v)
% where v is the perpendicular
%
% Example:

==> planar.m <==
% PLANAR: define a planar surface for ray tracing

==> plot_opt3d_surfaces.m <==
% PLOT_OPT3D_SURFACES: display surfaces used for 3d raytracing

==> plot_rays_for_length.m <==
% plot_rays_for_length: plot rays for a certain propagation distance

==> plot_rays_to_dest.m <==
% plot_rays_to_dest: plot rays as they head towards a particular point

==> points2ray_sinelens.m <==
% POINTS2RAY_SINELENS: solve for ray parameters given initial & final points

==> prepare_ideal_lens.m <==
% PREPARE_IDEAL_LENS: configure an ideal lens for raytracing

==> projective_combine.m <==
% projective_combine: consolidate sequential projective transforms

==> projective_conjugate.m <==
% Calculate the conjugate point for a projective transformation

==> ray.m <==
% RAY: a structure defining a ray

==> ray_convergence_error.m <==
% RAY_CONVERGENCE_ERROR: compute aberrations given a set of rays

==> ray_flat_intersection.m <==
% ray_flat_intersection: compute intersection of rays with flat surface

==> raytrace.m <==
% RAYTRACE: trace light rays through refractive surfaces

==> raytrace_ideal_lens.m <==
% RAYTRACE_IDEAL_LENS: ray propagation through well-corrected lenses

==> raytrace_matrix.m <==
% raytrace_matrix: a collection of utilities for matrix optics

==> raywaist.m <==
% RAYWAIST: calculate the width of a collection of rays at various

==> run_lightsheet_tracing.m <==
% unknown

==> Sellmeier.glass.refr <==
% database of the Sellmeier coefficients of glass meterials, useful to calculate 
% the refractive index for different wavelengths, i.e., chromatic
% dispersion

==> spherical.m <==
% unknown

==> symsph.m <==
% SYMSPH: symmetric spherical lens

==> trace_mirror_zernike.m <==
% TRACE_MIRROR_ZERNIKE: trace rays reflecting from a mirror whose shape is specified by zernike coefficients

==> trace_rays_sineobj.m <==
% TRACE_RAYS_SINEOBJ: ray-tracing through an objective lens that satisfies the sine condition.

==> trace_rays_sineobjtube.m <==
% TRACE_RAYS_SINEOBJTUBE: ray-tracing through a microscope that satisfies the sine condition.

==> transmitted_ray.m <==
% TRANSMITTED_RAY: calculate the direction of transmitted ray

==> tune_dm1.m <==
% unknown

==> tune_dm1_linear.m <==
% unknown

==> tune_dm2d.m <==
% unknown

==> tune_dm_bruteforce.m <==
% TUNE_DM: optimize the parameters of a deformable mirror

Subfolder Old
No any file

Subfolder unit_testing
==> test_ideal_lens.m <== 
% unknown