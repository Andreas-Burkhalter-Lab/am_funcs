function structout = pd_tune_single_actuator(im,im_DM,v,pupildata,options)
%   If you plan to do a Zernike-basis optimization (see "do_zernike"
%   below), then options must contain the field
%     Zindex: a sorted list of Zernike basis functions (e.g., 1:20 for all
%       up through fifth order)
%   Additionally, options may have the following fields:
%     do_gaussian (default true): if true, the user will have the chance to
%       define, visualize, and/or improve a gaussian approximation to the
%       aberration
%     do_zernike (default true): if true, the user will be asked to 
%       define, visualize, and/or improve a Zernike-basis approximation to the
%       aberration
%     gaussianP: starting guess for the gaussian aberration (see
%       pd_gaussian_gui for format)
%     phiV: an explicit representation of the phase from the Gaussian fit,
%       can be supplied as an initializer to the Zernike optimization if
%       you want to bypass the Gaussian step
%     Zvalue, Zxreg, Zreg: these "outputs" (see below) can be passed as
%       inputs for further optimization
%     interactive (default true): if true, runs GUIs in blocking mode
%       waiting for user input (if false, performs optimizations
%       automatically)

% Copyright 2009 by Timothy E. Holy

  if nargin < 5
    options = struct;
  end
  if isfield(options,'gaussianP') && isempty(options.gaussianP)
    options = rmfield(options,'gaussianP');
  end
  options = default(options,'do_gaussian',true,'do_zernike',true,'interactive',true,'gaussianP',[200 0 0 1/0.6^2]);
  if ~isequal(size(im),size(im_DM))
    error('The two image stacks must be of the same size');
  end
  if options.do_zernike
    if (length(intersect(options.Zindex, [1 2])) ~= 2)
      error('You must have both first-order Zernikes in Zindex');
    end
    if ~issorted(options.Zindex)
      error('Zindex must be sorted');
    end
  end
  structout = rmfield(options,{'do_gaussian','do_zernike','interactive'});
  if options.do_gaussian
    % Get an initial guess from the user for the voltage-dependent aberration
    ops.paramsin = options.gaussianP;
    ops.interactive = options.interactive;
    if options.interactive
      [structout.gaussianP,structout.phiV] = pd_gaussian_gui(im_DM,v,pupildata,ops);
    else
      hfig = pd_gaussian_gui(im_DM,v,pupildata,ops);
      handles = guidata(hfig);
      pd_gaussian_gui('btnOptimize_Callback',hfig, [], handles);
      [structout.gaussianP,structout.phiV] = pd_gaussian_gui('btnDone_Callback',hfig, [], handles);
      close(hfig);
    end
    if isempty(structout.phiV)
      % User closed the figure, cancel
      return
    end
    structout.gaussianGUI = true;  % a flag to indicate that the user used gaussianGUI
  end
  if options.do_zernike
    if ~isfield(structout,'Zvalue') || isempty(structout.Zvalue)
      % Project phiV to a Zernike basis to intialize Zvalue
      zernv = zernike_values(pupildata.rho,pupildata.theta,options.Zindex);
      N = sum(pupildata.H0(:));
      nZ = length(options.Zindex);
      Zv = zeros(1,nZ);
      for zi = 1:nZ
        tmp = structout.phiV .* zernv(:,:,zi);
        Zv(zi) = sum(tmp(:))/N;
      end
      % Construct a linear-model guess for the voltage-dependent aberration
      structout.Zvalue = v(:) * Zv;
    end
    % Pick the base image
    [tmp,baseIndex] = mindist(0,v);
    % Do a coarse cross-registration (so there is enough overlap for
    % fine-scale matching)
    if ~isfield(structout,'Zxreg')
      imbase = imfilter_gaussian(im(:,:,baseIndex),[3 3]);
      imDMbase = imfilter_gaussian(im_DM(:,:,baseIndex),[3 3]);
      structout.Zxreg = register_zernike({imbase,imDMbase},pupildata);
    end
    if isfield(structout,'Zreg')
      % Re-introduce the registration compensation
      structout.Zvalue(:,[1 2]) = structout.Zvalue(:,[1 2]) + structout.Zreg;
    end
    % Optimize these values
    ops.interactive = options.interactive;
    ops.Zxreg = structout.Zxreg;
    if options.interactive
      structout.Zvalue = pd_zernike_gui(im(:,:,baseIndex),im_DM,pupildata,options.Zindex, ...
        structout.Zvalue,ops);
    else
      hfig = pd_zernike_gui(im(:,:,baseIndex),im_DM,pupildata,options.Zindex, ...
        structout.Zvalue,ops);
      handles = guidata(hfig);
      fval = pd_zernike_gui('btnOptimize_Callback',hfig, [], handles);
      fval = fval(end);
      fvalOld = fval+1;
      while (fval < fvalOld)
        % Keep doing linear regression & re-fitting until it stops getting
        % better
        fvalOld = fval;
        pd_zernike_gui('btnLinRegress_Callback',hfig, [], handles);
        fval = pd_zernike_gui('btnOptimize_Callback',hfig, [], handles);
        fval = fval(end);
      end
      fvalOld = fval+1;
      while (fval < fvalOld)
        % Keep doing median filtering & re-fitting until it stops getting
        % better
        fvalOld = fval;
        pd_zernike_gui('btnMedian_Callback',hfig, [], handles);
        fval = pd_zernike_gui('btnOptimize_Callback',hfig, [], handles);
        fval = fval(end);
      end
      structout.Zvalue = pd_zernike_gui('btnDone_Callback',hfig, [], handles);
      close(hfig);
    end
    % Register the unaberrated stack to account for any drift components
    if ~isfield(structout,'Zreg')
      pupildata.baseIndex = baseIndex;
      structout.Zreg = register_zernike(im,pupildata);
    end
    if ~isempty(structout.Zvalue)
      structout.Zvalue(:,[1 2]) = structout.Zvalue(:,[1 2]) - structout.Zreg;
      structout.zernikeGUI = true;
    else
      % User closed without hitting "Done" (so it's a cancel)
      structout.zernikeGUI = false;
    end
  end
  