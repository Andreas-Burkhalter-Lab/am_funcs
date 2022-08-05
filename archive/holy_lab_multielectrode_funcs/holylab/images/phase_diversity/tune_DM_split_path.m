function pout = tune_DM_split_path(image_pairs,pair_voltages,pupil_data,p0,mirror_func,options)
% image_pairs: h-by-w-by-2-by-n_pairs (first of pair is
%   "unaberrated")
% pair_voltages: n_actuators-by-n_pairs set of mirror voltage parameters
% pupil_data:
% p0: a vector of parameters to mirror_func that describe other
%   components (typically the actuator geometry); these are the ones that
%   will be tuned
% mirror_func: a function handle that returns the mirror shape,
%      z = mirror_func(x,voltages,p)
% where shape is evaluated at positions x
  
% Note: this function isn't finished!
  
  %% Argument parsing
  impsz = size(image_pairs);
  if (length(impsz) ~= 4)
    error('Image data must be 4-dimensional')
  end
  if (impsz(3) ~= 2)
    error('Must be supplied as image pairs');
  end
  n_pairs = impsz(end);
  
  %% Do the pupil stuff here to get pupil coords
  % Load the pupil coords that matter (i.e., the ones inside the pupil)
  %   into the n_pts-by-2 matrix X
  % Create H0
  % Create pupil_index as the index of points in X that are within the
  %   pupil of H0 (like the example in the help of mirao52)
  
  % fminsearch isn't ideal, but there aren't that many parameters to
  % optimize, so it should be OK
  pout = fminsearch(mirror_param_error,p0);
  
  function err = mirror_param_error(p)
    % written as a nested function
    err = 0;
    for i = 1:n_pairs
      z_mirror = mirror_func(X,pair_voltages(:,i),p);
      phi = zeros([size(H0) 2]);
      phi(pupil_index+numel(H0)) = (2*pi/lambda) * z_mirror;
      err = err + pdpenalty(phi,image_pairs(:,:,:,i),H0);
    end
  end
end
