function [pout,imc] = pd_analyze_vdep_stack(im_ab,v,pupildata,params)
% PD_ANALYZE_VDEP_STACK: extract the voltage-dependent component of the aberration by phase diversity
% 
% Note: this is a GUI, so it requires user input. The GUI will play the
% acquired frames, do a preliminary pairwise fit, and then ask the user
% to define the range of frames to be used in the linear fit. The user
% specifies this range by sliding the vertical lines to the start and end
% of the range and then clicking "Continue".
%
% Syntax:
%   [p,imc] = pd_analyze_vdep_stack(im_ab,v,pupildata,params)
% where
%   im_ab is an m-by-n-by-K array, so that im_ab(:,:,k) is the kth image
%     in a sequence;
%   v is a vector of length K, v(k) is the voltage (or amplitude)
%     corresponding to the kth image
%   pupildata is a structure which must have the following fields:
%     H0, rho, theta: pupil parameters
%   and params is a structure which may have the following fields:
%     Registration-related:
%       Zreg: a K-by-2 matrix, listing the first-order Zernike coefficients
%         that register the unaberrated (i.e., DM-free) "channel"
%       OR
%       im_unab: m-by-n-by-K array, containing the image data from the
%         unaberrated channel
%       If neither of these is supplied, the algorithm assumes the images are
%       already registered (i.e., Zreg = 0).
%     Zvdep: the initial guess for the voltage-dependent aberration,
%       expressed as coefficients of the Zernike basis (functions indexed by
%       Zindex). If this is not supplied, an initial guess will be
%       generated.
%     Zindex (default 1:20): the list of Zernike polynomials (in
%       single-index notation) to be used to fit the images
%     selIndex: a vector specifying particular frames to be used in fitting
%       the voltage-dependent phase (default: will ask the user in a GUI
%       format)
%     use_raw (default false): if true, fits the final phase as a "raw" matrix
%       rather than using Zernike polynomials
%     baseIndex: the image # to be used as the "base" image (default: the
%       one closest to v = 0)
%     display (default true): if true, plots phase images on the screen
%     sort_v (default true): if true, will internally re-order the
%       voltages to be sorted (if not already sorted). One case in which
%       this might not be desirable is if you oscillate the voltage and
%       there is substantial drift
% and
%   p is a structure containing the following fields:
%     selIndex: a vector containing the frame numbers/actuator voltages
%       used to compute the voltage-dependent phase.
%     Zreg: the first-order Zernike coefficients used for registration
%       (will be zero unless params.Zreg or params.im_unab was provided)
%   and EITHER
%     phiV: the raw voltage-dependent phase
%   OR
%     Zindex,Zvalue: vectors listing the Zernike basis functions (in
%       single-index notation) and their coefficient values
%  imc is an m-by-n-by-K array containing the computed images for the
%    frames listed in p.selIndex.
%
% It's worth noting that the output "p" is already in a form suitable to
% serve as the foundation for a new "params" input. Thus, you can easily
% fit using a Zernike basis, and then transition to raw.
%
% See also: PD_ANALYZE_VDEP_STACK_DEMO.
  
% Copyright 2009 by Timothy E. Holy
  
  %% Argument parsing
  [mind,baseIndex] = mindist(0,v(:)');  % find frame with voltage closest to 0
  params = default(params,'sort_v',true,'Zindex',1:20,'baseIndex',baseIndex,'display',true,'use_raw',false);
  selIndex = 1:length(v);
  have_selIndex = isfield(params,'selIndex');
  if have_selIndex
    selIndex = params.selIndex;
    if isempty(intersect(params.baseIndex,selIndex))
      error('The base frame is not part of the selected set of frames')
    end
  end
  
%   if ~issorted(v) && params.sort_v
%     [v,sortOrder] = sort(v);
%     im_ab = im_ab(:,:,sortOrder);
%   else
%     sortOrder = 1:length(v);
%   end
  
  %% Handle registration, either as supplied or from DM-free images
  if isfield(params,'Zreg')
    Zreg = params.Zreg;
  elseif isfield(params,'im_unab')
    ptmp = params;
    % Do translation-based registration. Remove Zindex and Zc0 so that it
    % just uses the first-order Zernikes
    if isfield(ptmp,'Zindex')
      ptmp = rmfield(ptmp,'Zindex');
    end
    if isfield(ptmp,'Zc0')
      ptmp = rmfield(ptmp,'Zc0');
    end
    Zreg = register_zernike(params.im_unab,pupildata,ptmp); % note we deliberately register all frames, not just selIndexed ones
  else
    Zreg = zeros(length(v),2);
  end
  pout = struct('Zreg',Zreg);

  
  %% Create the initial guess for the voltage-dependent component of the aberration
  % If an initial guess for Zvdep was supplied, use that
  if isfield(params,'Zvdep')
    Zvdep = params.Zvdep;
  else
    % We have to calculate the starting guess for Zvdep
    % Do pairwise comparisons to the base frame
    if (have_selIndex && ~isequal(selIndex,1:length(v)))
      % if we're here, we are surely doing all of them
      warning('pd:ignore','Ignoring the user-supplied selIndex, because there was no Zvdep');
    end
    Zc = register_zernike(im_ab,pupildata,params);
    Zc(:,1:2) = Zc(:,1:2) - Zreg;

    % Ask the user to define the range (we don't do both ranges in this
    % case)
    hfig = figure;
    plot(v,Zc)
    yl = get(gca,'YLim');
    hl1 = line(v([2 2]),yl,'Color','r');
    hl2 = line(v([end-1 end-1]),yl,'Color','r');
    drag_line(hl1);
    drag_line(hl2);
    uicontrol('Parent',hfig,'Position',[20 20 200 40],'String','Continue','Callback','uiresume(gcbf)');
    uiwait(hfig);
    xd1 = get(hl1,'XData');
    xd2 = get(hl2,'XData');
    xd = [xd1(:); xd2(:)];
    close(hfig);
    selFlag = v >= min(xd) & v <= max(xd);
    selIndex = find(selFlag);

    % Do linear interpolation over the defined range
    n_terms = size(Zc,2);
    Zvdep = zeros(1,n_terms);
    Zstatic = zeros(1,n_terms);
    for k = 1:n_terms
      [Zvdep(k),Zstatic(k)] = linregress(v(selFlag),Zc(selFlag,k));
    end
  end
  pout.selIndex = selIndex;

  %% Compute the (starting) phase from the voltage-dependent Zernike basis
  zernv = zernike_values(pupildata.rho,pupildata.theta,params.Zindex);
  phiV0 = zeros(size(zernv(:,:,1)));
  for k = 1:length(params.Zindex)
    phiV0 = phiV0 + zernv(:,:,k) * Zvdep(k);
  end
  if params.display
    figure; imagesc(fftshift(phiV0)); title('phiV0'); colorbar; drawnow
  end
  
  %% Optimize the voltage-dependent phase using all selected images simultaneously
  % Set the mode and initial guess
  if params.use_raw
    params.mode = 'registered vdep raw singlechannel';
    p0 = phiV0;
  else
    params.mode = 'registered vdep Zernike singlechannel';
    p0 = Zvdep;
  end
  % Optimize the phase
  params.v = v(selIndex);
  params.Zshifts = Zreg(selIndex,:);
  zf = phi_parametrizations(params,pupildata);
  if params.use_raw
    % Limit changes to 1 (on a phase scale of 2*pi) in each iteration
    zf.normalize_grad = true;
    zf.mu_max = 1/max(abs(params.v));
  end
  zf.iter_max = 1000;
  [p,pout.fval] = calcphi2d(im_ab(:,:,selIndex),pupildata.H0,p0,zf);
  [p,tmp] = calcphi2d(im_ab(:,:,selIndex),pupildata.H0,p,zf); % a second to protect against false exits
  pout.fval = [pout.fval tmp];
  
  % Generate the calculated images
  if (nargout > 1)
    phik = zf.param2phi(p);
    [val,grad,object] = pdpenalty(phik,im_ab(:,:,selIndex),pupildata.H0);
    [Hk,sk,imc] = pd_forward_model_2d(phik,pupildata.H0,object);
%       selIndex{groupIterator} = sortOrder(selFlag{groupIterator});
  end
  
  % Display and munge output into proper format
  if params.use_raw
    phiV = p;
    pout.phiV = p;
  else
    phiV = zeros(size(zernv(:,:,1)));
    for k = 1:length(params.Zindex)
      phiV = phiV + zernv(:,:,k) * p(k);
    end
    pout.Zindex = params.Zindex;
    pout.Zvalue = p;
  end
  if params.display
    figure; imagesc(fftshift(phiV)); title('phiV'); colorbar; drawnow
  end
end
