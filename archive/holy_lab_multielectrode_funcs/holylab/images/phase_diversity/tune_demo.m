% The most important step is to provide a reasonable starting guess for
% A, sigma, and v2z (see help in mirao52). The idea would be that
% starting values would come from an analysis of a couple of
% pairs. Alternatively, since you can almost "see" where the actuators
% are in the defocused images, you can probably make some pretty good
% guesses without having to do much in the way of analysis (although keep
% in mind that you need to map the coordinates into the fourier plane,
% not the image plane---by defocusing you start to essentially image the
% fourier plane, so the correspondence is "simple" but there will be an
% overall scaling & translation that will be important to get in the
% right ballpark).

s0 = struct('A',[starting guess here],'sigma',scalar_value,'v2z',scalar_value,'v0',0);
[p0,fields,field_shape,sbase] = extract_fields(s0,'v2z','sigma','A');
% I did them in the order v2z, sigma, A because I think v2z is going to
% be the hardest to guess well on, but since it's a scalar it can be
% optimized easily. I think fminsearch starts with the first variable,
% then the second, etc, so I think these two scalars will be pretty
% well-set by the time it gets to A. Conversely, if they're way off, I
% could imagine it might start "messing up" a fairly good initial guess
% on A.

fill_func = @(p) fill_fields(fields,field_shape,p,sbase);
mirror_func = @(X,v,p) mirao52(X,v,fill_func(p));

% Now load all your image data, specify pupil parameters, etc
% Note tune_DM_split_path has been designed to use your in-focus images
% It will take more thinking to do the out-of-focus ones---even though
% those look very informative---and I think that analysis will be more
% sensitive to the z-position of the beads so might require a full 3d
% model. In contrast, with the in-focus images, you're right at the
% optimium of the focus, and so z-position only comes in quadratically
% and for small z-deviations can hopefully be ignored.
p = tune_DM_split_path(image_pairs,pair_voltages,pupil_data,p0, ...
		       mirror_func);

mirror_structure = fill_func(p);
% That's it! We'd save mirror_structure and use it until you shift or
% change the hardware.
% One might eventually want to tune v0 for each actuator
% (to "truly flatten" the mirror), but that won't be as hard, or as
% important, as getting these position & scaling parameters nailed down.
%
% But, at this point we're basically ready to program in
% particular aberrations and see how it changes the images.
