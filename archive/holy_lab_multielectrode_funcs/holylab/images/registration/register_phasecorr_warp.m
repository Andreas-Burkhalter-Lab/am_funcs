function imw = register_phasecorr_warp(u,im,options)
% register_phasecorr_warp: deform an image
% 
% Syntax:
%   imw = register_phasecorr_warp(u,im,options)
% where
%   u is the deformation (see register_phasecorr_improve or
%     register_phasecorr_refine);
%   im is the input image (see register_phasecorr_initialize for valid
%     choices) 
%   options comes from register_phasecorr_initialize. You can add one extra
%     field, 'shift', to supply an additional translation to the entire
%     image.  The units of the shift should be pixels, the same as for u.
% and
%   imw is the warped image.
%
% See also: register_phasecorr_initialize, register_phasecorr_improve, register_phasecorr_refine.
  
% Copyright 2010 by Timothy E. Holy
  
  uh = register_phasecorr_prolong_fullsize(u,options);
  if isempty(uh)
    imw = im;
    return
  end

  %% Calculate the warped image
  g = cell(1,options.n_dims);
  for dimIndex = 1:options.n_dims
    g{dimIndex} = options.g0{dimIndex} + uh{dimIndex};
  end
  imw = zeros(size(im),class(im));
  colons = repmat({':'},1,options.n_dims);
  for chanIndex = 1:options.n_values
    imw(colons{:},chanIndex) = iminterp(im(colons{:},chanIndex),g{:});
  end
  