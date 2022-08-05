% There are several distinct frameworks for registration.
%
% ::Top-level wrappers/GUIs::
% These functions may be useful for processing a series of stacks (image
% volumes).
% They call elements of the core modules described below.
% multigrid_registration_stepper  - multigrid_vcycle-based wrapper for various registration methods
% jpm_translation_registration_script - Registering stacks
% register_dfof_validate_stacks   - a GUI for checking for bad frames in base stacks
% register_movie_dfof_gui         - register a movie using stimulus responses
% register_movie_gui              - register grayscale volumes over time
% register_movie_smoothu_gui      - enforce deformation smoothness over time by filtering
% register_options_dialog         - choosing options for GUI registration
%
%
% ::Core modules::
%
% Translational-only (rigid) registration:
% register_translate_nancorr - align images by minimizing mismatch (recommended)
% register_translate_pad     - pads an image with NaNs in preparation for nancorr
% register_translate         - align images by phase correlation (alternative to nancorr)
% register_rigid             - rigid registration (limited functionality alternative, not recommended)
%
%
% Rigid registration using control points:
% cp2tform_3d_affine         - This func return the affine tform
% cp2tform_3dsimilarity      - compute a three-dimensional similarity transform from control points
%
% Warp registration (non-rigid) by multigrid:
% register_demo	              - This is a short demo script on running multigrid registration.  There is no help for this, you simply need to look at the code.
% register_multigrid_vcycle   - non-rigid image registration using a multigrid algorithm
% register_multigrid_options  - initialize multigrid registration
% register_change_level_gap   - change scale of deformation grid in registration
% register_gradu_restrict     - move gradients with respect to u to a coarser grid
% register_mismatch_noncov    - compute the mismatch error under deformation
% register_multigrid_penalty  - penalty, gradient, and hessian for multigrid registration
% register_multigrid_auto     - non-rigid image registration using a multigrid algorithm
% register_rigid2nonrigid     - convert translation into a multigrid-able deformation
% register_shift2mask         - use parameters of translation to estimate a mask
% stack_automask              - takes an input img and creates a custom 3d mask
% drawmask_gui                - specify mask interactively
% register_multigrid_options_compress - strip out big entries for disk-storage
% register_multigrid_options_expand   - restore the full parameter structure
% uinterp                     - interpolates "u" values from image registration (multi-stack wrapper)
% warpnwrite                  - warpnwrite will warp and write a .cam file to an already-open file
% fill_bad_frames             - attempt to identify and interpolate bad frames from a .imagine/.cam field
% padnans                     - pads an image with the edges supplied in nanpad
% 
%
% Warp registration by block matching:
% register_block_improve    - discover or refine a deformation to improve image registration
% register_block_penalty    - compute the data (mismatch) + smoothness penalty function
% register_block_mismatch   - compute image mismatch in tapered blocks
% register_mismatch_nancorr - compute mismatch under translation, allowing missing data (similar to register_translate_nancorr)
%
%
% Warp registration by phase correlation:
% register_phasecorr_composeu 		- compose deformations for phasecorr registration
% register_phasecorr_improve		- optimize a deformation for registering images
% register_phasecorr_initialize		- set parameters for phasecorr registration
% register_phasecorr_multires		- iterative refinement of deformation until improvement stops
% register_phasecorr_prolong		- resize to a finer grid
% register_phasecorr_prolong_fullsize - bring u to full size of image
% register_phasecorr_refine			- change a coarse deformation to a finer one
% register_phasecorr_warp			- deform an image
%
%
%
%
% ::Utility functions::
% These are used in more than one framework.
% register_movie_apply    - apply pre-calculated deformations to a sequence of stacks, to make a registered sequence
% register_multigrid_warp - convenience function for nonrigid image deformations
% register_composeg       - Compose one deformation with another
% register_g0             - create default deformation (identity map)
% g_array2cell            - convert array-based representation into cell array
% register_u_prolong      - move u to a finer grid
% register_u_restrict     - transfer u to a coarser grid
% register_visualize_u    - view a deformation in 2d or 3d
% nanbox                  - find the largest subarray containing no NaNs
% register_logdetpenalty  - penalize deformations based on log(det(Jacobian))
% register_logdetpenalty_composition - explicit composition with a log(det(J)) penalty
% register_gui_utilities  - a collection of function handles often used in registration, particularly by the GUI programs.
% 


