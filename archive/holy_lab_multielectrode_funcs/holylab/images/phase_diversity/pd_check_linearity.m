function [v,Z,fval_all] = pd_check_linearity(imdata,pupildata,act_fit)
% PD_CHECK_LINEARITY: fit Zernikes to individual frames
% Syntax:
%   [v,Z,fval] = pd_check_linearity(imdata,pupildata,act_fit)
% where act_fit contains the actuator fitting data for a single actuator
% assuming a piecewise linear model (e.g., as output from
% PD_TUNE_MIRROR_FROM_GRID).
%
% See also: PD_TUNE_MIRROR_FROM_GRID.
  
% Copyright 2009 by Timothy E. Holy
  
  actuator = unique([act_fit.actuator]);
  if (length(actuator) ~= 1)
    error('Must supply fit information from a single actuator');
  end
  actuatorFlag = imdata.actuator == actuator;
  K = length(imdata.v);
  n_fits = length(act_fit);
  [tmp,baseIndex] = mindist(0,imdata.v);
  v = [];
  nZ = length(act_fit(1).Zindex);
  Z = zeros(0,nZ);
  fval_all = zeros(0,n_fits);
  for k = 1:K
    if (k == baseIndex)
      continue;
    end
    if isfield(pupildata,'snipindx')
      imtmp = imdata.stack(pupildata.snipindx{:},[baseIndex k], ...
        actuatorFlag);
    else
      imtmp = imdata.stack(:,:,[baseIndex k],actuatorFlag);
    end
    thisv = imdata.v(k);
    fval = zeros(1,n_fits);
    Zfit = zeros(n_fits,nZ);
    for fitIndex = 1:n_fits
      Zc0 = act_fit(fitIndex).Zvalue * thisv;
      [Zfit(fitIndex,:),fval(fitIndex)] = register_zernike(imtmp,pupildata,...
        struct('Zindex',act_fit(fitIndex).Zindex,'Zc0',Zc0));
    end
    [fvaltmp,bestFitIndex] = min(fval);
    v(end+1) = thisv; %#ok<AGROW>
    Z(end+1,:) = Zfit(bestFitIndex,:); %#ok<AGROW>
    fval_all(end+1,:) = fval; %#ok<AGROW>
  end
end
