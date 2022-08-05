function mask = register_shift2mask(sz,shift,options)
  if (nargin < 3)
    options = struct;
  end
  options = default(options,'margin_frac', 0.1, 'margin_pixel', 2);
  
  sz = sz(1:size(shift,2));
  mask = ones(sz);
  amplitude = 1+options.margin_frac;
  minshift = min(shift,[],1);
  mask_shift = image_shift(mask, floor(amplitude*minshift - options.margin_pixel));
  mask(isnan(mask_shift)) = 0;
  maxshift = max(shift,[],1);
  mask_shift = image_shift(mask, ceil(amplitude*maxshift + options.margin_pixel));
  mask(isnan(mask_shift)) = 0;
end