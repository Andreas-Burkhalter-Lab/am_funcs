function [out1,out2] = ocpi_pdwrapper(varargin)
% OCPI_PDWRAPPER: homogenous phase-diverse OCPI.
% This can handle a gaussian-profile illumination beams (assumed
% constant across the field of view) and aberrations expressed in terms
% of Zernike coefficients, as well as arbitrary functional forms.
% Syntax:
%   s = ocpi_pdwrapper(inputparams,outputparams)
% An all-structures syntax; most useful if you want this to do some
% calculations for you.
%   [val,grad] = ocpi_pdwrapper(p,inputparams)
% Used if you want to make this the penalty function for a minimization.
%   s = ocpi_pdwrapper(p,inputparams,'insert')
% Do this to re-sync the parameters after an optimization.
%   p = ocpi_pdwrapper(inputparams,'extract')
%   p = ocpi_pdwrapper(inputparams,'extract',fieldnames)
% Do this to extract a coordinate vector for use by optimization routines.
%
% inputparams fields:
%   rho: pupil data
%   theta: "
%   H0:    "
%   c0: defocus pupil data
%   sigma OR v: parameters of illumination
%   phik OR Zinfo OR DMinfo: parameters of aberrations
%     Zinfo fields:
%       Zindex: a cell array of length K, each entry a vector of Zernike
%         coefficient indices for each diversity image
%       Zcoefs: a cell array of the same configuration as Zindex, giving
%         the coefficient values
%       map: a cell array of the same configuration as Zindex, where
%         each element is a vector containing the index of a parameter
%         vector used for optimization.  0 indicates that the parameter is
%         fixed (not optimized). This notation allows you to "tie"
%         particular coefficients together across diversity images (or you
%         can make them all separate, it's up to you).
%       Example:
%         Zindex = {7, [4 7 12]}
%         Zcoefs = {0.5, [10 0.5 1.3]}
%         map = {1, [0 1 2]}
%       would indicate that both diversity images have a common coma
%       aberration of 0.5 that is intended for optimization. The second
%       diversity image also has a fixed defocus of 10 and an optimizable
%       spherical aberration of 1.3.  The parameter p for optimization
%       would consist of 2 values (for the Zernikes, that is), where the
%       first parameter is the common value of coma and the second is the
%       spherical aberration applied to the 2nd diversity image.
%     DMinfo fields: TBD
%   f (the object) OR psi (sqrt of the object) (use the latter to force
%     nonnegative objects)
%   Iobsk (optional): the observed diversity stacks
%   mask: 0s or 1s depending on which pixels to keep
%   stacklen (needed if supplying neither psi, Iobsk, nor v): # of frames in
%     a stack
%
% If you're using the p,inputparams syntax, there's an extra field:
%   datafields: a cell array of field names, describing the order in
%     which the values in p should be placed into inputparams.
%
% outputparams fields (you probably don't have to mess with these
% directly, there is an attempt at intelligent defaults):
%   PSFs: if true, returns the PSFs for each diversity image
%   Iobs_calc: if true, returns the calculated images (the PSFs convolved
%     with the underlying object)
%   fgrad, vgrad, phikgrad: each determines whether gradients with
%     respect to the given parameter are calculated.
  

% An issue worth noting: this syntax probably isn't ideal for handling
% the case in which all K diversity images have a common "overall"
% aberration with only a few Zernikes that differ (e.g., as would be used
% in conventional PD).  There will surely need to be some additional work
% done for this case.  However, this current syntax should be very useful
% for calibrating the DM.
  
  persistent Zval ZvalIndex
  
  % Parse the input and set up the defaults
  outputparams = struct;
  insert = false;
  extract = false;
  if isstruct(varargin{1})
    input_p = false;
    inputparams = varargin{1};
    if (length(varargin) > 1)
      if isstruct(varargin{2})
        outputparams = varargin{2};
      elseif ischar(varargin{2})
        extract = strcmp(varargin{2},'extract');
      end
    end
  else
    input_p = true;
    p = varargin{1};
    inputparams = varargin{2};
    if (length(varargin) > 2)
      if ischar(varargin{3}) && strcmp(varargin{3},'insert')
        insert = true;
      end
    end
  end
  
  % Extracting parameters is easy, do that first
  if extract
    if (length(varargin) > 2)
      datafields = varargin{3};
    else
      datafields = inputparams.datafields;
    end
    p = zeros(0,1);
    for indx = 1:length(datafields)
      thisfield = datafields{indx};
      switch thisfield
        case {'f','psi','sigma','v','phik'}
          tmp = inputparams.(thisfield);
        case 'Zcoefs'
          Zc = cat(2,inputparams.Zinfo.Zcoefs{:});
          map = cat(2,inputparams.Zinfo.map{:});
          tmp = nan(1,max(map)+1);
          tmp(map+1) = Zc; % +1 to handle the 0 index
          tmp(1) = []; % to clear the 0 index
      end
      p(end+1:end+numel(tmp)) = tmp(:);
    end
    out1 = p(:);
    return
  end


  %compute_gradients_default = input_p && nargout > 1;
  compute_gradients_default = false;  % this is handled better below
  outputparams = default(outputparams,'fgrad',compute_gradients_default,...
    'vgrad',compute_gradients_default,...
    'phikgrad',compute_gradients_default);

  % Determine whether the underlying object is supplied, and if so in
  % what form
  have_object = false;
  if isfield(inputparams,'f')
    obj_mode = 'f';
    have_object = true;
  elseif isfield(inputparams,'psi')
    obj_mode = 'psi';
    have_object = true;
  end

  % If the user is supplying only a limited number of fields, set up to
  % return those things that can be calculated
  if ~have_object
    % We can only return the PSFs
    outputparams = default(outputparams,'PSFs',true);
  end
  if have_object
    if ~isfield(inputparams,'Iobsk')
      % We can't return the error, so return the calculated images
      outputparams = default(outputparams,'Iobs_calc',true);
    else
      outputparams = default(outputparams,'Iobs_calc',false);
    end
  end

  % Determine the sizes of objects: imsz & K
  % Also determine the image class, imclass
  K = [];
  if isfield(inputparams,'Iobsk')
    sz = size(inputparams.Iobsk);
    K = sz(4);
    imsz = sz(1:3);
    imclass = class(inputparams.Iobsk);
  elseif have_object
    imsz = size(inputparams.(obj_mode));
    imclass = class(inputparams.(obj_mode));
  else
    if isfield(inputparams,'v')
      stacklen = length(inputparams.v);
    else
      stacklen = inputparams.stacklen;
    end
    imsz = [size(rho) stacklen];
    imclass = class(rho);
  end
  if isempty(K)
    if isfield(inputparams,'phik')
      K = size(inputparams.phik,3);
    elseif isfield(inputparams,'Zinfo')
      K = length(inputparams.Zinfo.Zindex);
    end
  end
  if length(imsz) < 3
    imsz(3) = 1;
  end
  stacklen = imsz(3);
  framesz = imsz(1:2);
  % Ideally one would do some cross-checking here...

  % If needed, fill in the values of inputparams from p
  if input_p
    offset = 0;
    for indx = 1:length(inputparams.datafields)
      switch inputparams.datafields{indx}
        case {'f','psi'}
          npix = prod(imsz);
          tmp = p((1:npix) + offset);
          offset = offset + npix;
          tmp = reshape(tmp,imsz);
          inputparams.(inputparams.datafields{indx}) = tmp;
          outputparams.fgrad = nargout > 1;
        case 'Zcoefs'
          mapc = inputparams.Zinfo.map;
          map = cat(2,mapc{:});
          n = max(map);
          tmp = p((1:n)+offset);
          for indx2 = 1:K
            thismap = mapc{indx2};
            thismapI = thismap(thismap > 0);
            inputparams.Zinfo.Zcoefs{indx2}(thismap>0) = ...
              reshape(tmp(thismapI),size(thismapI));
          end
          offset = offset+n;
          outputparams.phikgrad = nargout > 1;
        case 'sigma'
          inputparams.sigma = p(1+offset);
          offset = offset+1;
          outputparams.vgrad = nargout > 1;
        case 'v'
          inputparams.v = p((1:stacklen)+offset);
          offset = offset + stacklen;
          outputparams.vgrad = nargout > 1;
        case 'phik'
          phiksz = [framesz K];
          npix = prod(phiksz);
          phik = p((1:npix) + offset);
          offset = offset + npix;
          inputparams.phik = reshape(phik,phiksz);
          outputparams.phikgrad = nargout > 1;
      end
    end
  end

  if insert
    out1 = inputparams;
    return
  end
  
  % If needed compute f from psi
  if strcmp(obj_mode,'psi')
    inputparams.f = inputparams.psi.^2;
  end

  % Determine the illumination mode, and if needed compute v
  z = 0:stacklen-1;
  zwrap = (z > stacklen/2);
  z(zwrap) = z(zwrap) - stacklen;
  if isfield(inputparams,'sigma')
    ill_mode = 'sigma';
    v = exp(-z.^2/(2*inputparams.sigma^2));
    v = v / sum(v);  % normalize
    inputparams.v = v;
  else
    % The user supplied v
    ill_mode = 'v';
  end

  % Determine the aberration mode, and if needed compute phik
  if isfield(inputparams,'phik')
    ab_mode = 'phik';
  else
    ab_mode = 'Zinfo';
    uZindex = unique(cat(2,inputparams.Zinfo.Zindex{:}));
    if (isempty(ZvalIndex) || ~isequal(uZindex,ZvalIndex) || ~isequal(size(Zval(:,:,1)),size(inputparams.rho)))
      % Initialize our computation of the Zernike polynomials
      ZvalIndex = uZindex;
      Zval = zernike_values(inputparams.rho,inputparams.theta,ZvalIndex);
      Zval = cast(Zval,imclass);
    end
    phiksz = [framesz K];
    phik = zeros(phiksz,imclass);
    for indx = 1:K
      Zindex = findainb(inputparams.Zinfo.Zindex{indx},ZvalIndex);
      phik(:,:,indx) = Zcoefs2phi(inputparams.Zinfo.Zcoefs{indx},...
        Zval(:,:,Zindex));
    end
    inputparams.phik = phik;
  end

  % Do the big calculation
  s = phasedivocpi_val_grad(inputparams,outputparams);

  % Convert the illumination gradient, if needed
  if (strcmp(ill_mode,'sigma') && outputparams.vgrad)
    zvar = sum(z.^2 .* inputparams.v);
    dvdsigma = (z.^2 - zvar) .* inputparams.v / inputparams.sigma^3;
    s.sigmagrad = sum(s.vgrad .* dvdsigma);
  end

  % Convert the aberration gradient, if needed
  if (strcmp(ab_mode,'Zinfo') && outputparams.phikgrad)
    map = cat(2,inputparams.Zinfo.map{:});
    Zindex = cat(2,inputparams.Zinfo.Zindex{:});
    Kmap = cell(1,K);
    for indx = 1:K
      Kmap{indx} = repmat(indx,[1 length(inputparams.Zinfo.Zindex{indx})]);
    end
    Kmap = cat(2,Kmap{:});
    clabel = agglabel(map+1); % +1 to handle the 0 map
    clabel(1) = [];
    n = max(map);
    s.Zgrad = zeros(n,1);
    for indx = 1:n
      thisindx = clabel{indx};
      for indx2 = 1:length(thisindx)
        thisindx2 = thisindx(indx2);
        % Project onto Zernike polynomials
        tmp = s.phikgrad(:,:,Kmap(thisindx2)) .* Zval(:,:,find(ZvalIndex == Zindex(thisindx2)));
        s.Zgrad(indx) = s.Zgrad(indx) + sum(tmp(:));
      end
    end
  end

  % If needed, pack the value and gradients into the 2-parameter return
  if input_p
    out1 = s.val;
    if (nargout > 1)
      out2 = p;
      offset = 0;
      for indx = 1:length(inputparams.datafields)
        switch inputparams.datafields{indx}
          case {'f','psi'}
            npix = prod(imsz);
            tmp = s.fgrad(:);
            if strcmp(obj_mode,'psi')
              % Use chain rule to handle sqrt(object) input
              tmp = 2*inputparams.psi(:) .* tmp;
            end
            out2((1:npix) + offset) = tmp;
            offset = offset + npix;
          case 'Zcoefs'
            n = length(s.Zgrad);
            out2((1:n)+offset) = s.Zgrad;
            offset = offset+n;
          case 'sigma'
            out2(1+offset) = s.sigmagrad;
            offset = offset+1;
          case 'v'
            out2((1:stacklen)+offset) = s.vgrad;
            offset = offset + stacklen;
          case 'phik'
            phiksz = [framesz K];
            npix = prod(phiksz);
            out2((1:npix) + offset) = s.phikgrad(:);
            offset = offset + npix;
        end
      end
    end
  else
    out1 = s;
    if (strcmp(obj_mode,'psi') && outputparams.fgrad)
      out1.psigrad = 2 * inputparams.psi .* out1.fgrad;
      out1 = rmfield(out1,'fgrad');
    end
  end
end


function phi = Zcoefs2phi(Zcoefs,Zval)
  n_coefs = length(Zcoefs);
  sz = size(Zval);
  phi = zeros(sz(1:2),class(Zval));
  for indx = 1:n_coefs
    phi = phi + Zcoefs(indx)*Zval(:,:,indx);
  end
end
