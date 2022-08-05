function imOut = padnans(imIn, nanpad)
% function padnans pads an image with the edges supplied in nanpad
% See also multigrid_registration_stepper
tmpimBase = cat(2,nanpad{1},imIn);
tmpimBase = cat(2,tmpimBase,nanpad{1});
tmpimBase = cat(1,nanpad{2},tmpimBase);
tmpimBase = cat(1,tmpimBase,nanpad{2});
tmpimBase = cat(3,nanpad{3},tmpimBase);
tmpimBase = cat(3,tmpimBase,nanpad{3});
imOut = tmpimBase;
end