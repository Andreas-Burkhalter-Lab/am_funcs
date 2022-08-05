function act_fit = pd_tune_mirror_from_grid(imdata,act_fit0,pupildata,act_geometry)
  % imdata: structure with fields
  %   stack: m-by-n-by-length(v)-by-n_actuators array of images
  %   v: voltages applied to individual actuators
  %   actuator: a list of actuator numbers
  % actdata: a structure array of the form output by pd_analyze_vdep_stack,
  %   with the additional field "actuator" telling which actuator is being
  %   used for each element. There should be only one per actuator. Note
  %   that the output of this function (act_fit) is also suitable for this
  %   input, as long as you pick one actuator.
  % pupildata: H0, rho, theta. May also have "snipindx", a 1-by-2 cell
  %   array containing sniprange with respect to imdata.
  % act_geometry: 
  %   amplitude (a scalar, try 100)
  %   centers: 52-by-2 array giving the coordinate centers (in pupil
  %     coordinates). See mirao52_layout, and transform with the tform
  %     created by pd_fitphase2grid_gui.
  %   sigma: scalar, the width of a given actuator. (try 0.5)
  
  actuators = [act_fit0.actuator];
  n_actuators = length(actuators);
  
  isInPupil = pupildata.rho <= 1;
  rho = pupildata.rho(isInPupil(:));
  theta = pupildata.theta(isInPupil(:));
  X = [rho .* cos(theta), rho .* sin(theta)];
  p2f = chainrule_gaussian(X,isInPupil,false);
  nH = sum(pupildata.H0(:));
  
  Zindex = unique([act_fit0.Zindex]);
  zernv = zernike_values(pupildata.rho,pupildata.theta,Zindex);
  
  % Do the fit twice, once for v <= 0 & once for v >= 0
  selIndex = {find(imdata.v <= 0),find(imdata.v >= 0)};
  
  act_fit = [];
  for actIterator = 1:n_actuators
    % Generate a gaussian-shaped initial guess
    thisActuator = actuators(actIterator);
    thisZindex = act_fit0(actIterator).Zindex;
    thisZvalue = act_fit0(actIterator).Zvalue;
    xy = act_geometry.centers(thisActuator,:);
    phiV = p2f([act_geometry.amplitude; xy(:); 1/act_geometry.sigma^2]);
    % Project onto a Zernike basis
    Zindex_index = findainb(thisZindex,Zindex);
    Zc0 = zeros(1,length(Zindex_index));
    for k = 1:length(Zindex_index)
      tmp = phiV .* zernv(:,:,Zindex_index(k));
      Zc0(k) = sum(tmp(:));
    end
    Zc0 = Zc0 / nH;
    % Compare these against the previous Zernikes (flipped into the proper
    % orientation)
    phiV0 = zeros(size(phiV));
    for k = 1:length(thisZvalue)
      phiV0 = phiV0 + thisZvalue(k) * zernv(:,:,Zindex_index(k));
    end
    tmp = phiV .* phiV0;
    s = sum(tmp(:));
    phiV0sym = fliplr(phiV0); phiV0sym = -flipud(phiV0sym);
    tmp = phiV .* phiV0sym;
    if (sum(tmp(:)) > s)
      % Flip it
      for k = 1:length(thisZvalue)
        if (mod(zernfun2(thisZindex(k)),2) == 0)
          thisZvalue(k) = -thisZvalue(k);
        end
      end
    end
    % Replace the first-order Zernike coefficients with those from the
    % previous fit (the mirror has more tip/tilt than predicted by the
    % gaussian)
    firstOrderFlag = thisZindex < 3;
    Zc0(firstOrderFlag) = thisZvalue(firstOrderFlag);
    
    % Do the fit(s)
    stackIndex = find([imdata.actuator] == thisActuator);
    if isfield(pupildata,'snipindx')
      im_ab = imdata.stack(pupildata.snipindx{:},:,stackIndex);
    else
      im_ab = imdata.stack(:,:,:,stackIndex);
    end
    fit0 = act_fit0(actIterator);
    fit0.Zvdep = Zc0;
    fit0.display = false;
    n_fits = length(selIndex);
    % First we fit with the gaussian starting guess for each subset of
    % frames
    for k = n_fits:-1:1
      fit0.selIndex = selIndex{k};
      results(k,k) = pd_analyze_vdep_stack(im_ab,imdata.v,pupildata,fit0);
    end
    % Then we fit each subset using the results for another subset
    for k = 1:n_fits
      for k1 = 1:n_fits
        if (k ~= k1)
          fit0.selIndex = selIndex{k};
          fit0.Zvdep = results(k1,k1).Zvalue;
          results(k,k1) = pd_analyze_vdep_stack(im_ab,imdata.v,pupildata,fit0);
        end
      end
    end
    % Finally, fit each subset using the old results
    for k = 1:n_fits
      fit0.selIndex = selIndex{k};
      fit0.Zvdep = thisZvalue;
      results(k,n_fits+1) = pd_analyze_vdep_stack(im_ab,imdata.v,pupildata,fit0);
    end
    % Now choose the best results for each subset
    fval = zeros(1,n_fits);
    for k = 1:n_fits
      for k1 = 1:n_fits+1
        fval(k1) = results(k,k1).fval(end);
      end
      [minval,minIndex] = min(fval);
      new_results = results(k,minIndex);
      new_results.actuator = thisActuator;
      if isempty(act_fit)
        act_fit = new_results;
      else
        act_fit(end+1) = new_results;
      end
    end
  end
end