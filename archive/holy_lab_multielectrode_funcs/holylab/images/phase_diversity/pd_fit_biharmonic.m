function params = pd_fit_biharmonic(data_specs,imDM,internalData,params0,Zreg,options)
  options = default(options,'debug',false,'Display','iter');
  if options.debug
    options = default(options,'MaxFunEvals',1);
  end
  % params = pd_fit_biharmonic(data_specs,camDMfilename,snipindx,internalData,params0)
  fn = fieldnames(params0);
  [p0,fields,field_shape,sbase] = extract_fields(params0,fn{:});
  fillfunc = @(p) fill_fields(fields,field_shape,p,sbase);

  if (size(Zreg,2) == 3)
    Zdefocus = Zreg(:,3);
  elseif (size(Zreg,2) == 1)
    Zdefocus = Zreg;
  else
    error('Zreg format not recognized');
  end

  optfunc = @(p) pd_err_biharm_total(p,imDM,data_specs,fillfunc,internalData,Zdefocus,options);
  p = fminsearch(optfunc,p0,options);
  
  params = fillfunc(p);
end
  
function err = pd_err_biharm_total(p,imDM,data_specs,fillfunc,internalData,Zdefocus,options)
  params = fillfunc(p);
  err = 0;
  n_act = size(internalData.act_xy,1);
  H0 = double(internalData.inpupil);
  for indx = 1:length(data_specs.act_num)
    act_num = data_specs.act_num(indx);
    v = zeros(1,n_act);
    K = length(data_specs.vol_tune);
    phi = zeros([size(internalData.inpupil) K]);
    for k = 1:K
      v(act_num) = data_specs.vol_tune(k);
      phi(:,:,k) = biharmonic_phi(params,v,internalData);
    end
    % Adjust the defocus
    Z4 = zeros(size(H0));
    Z4(internalData.inpupil) = internalData.Z4;
    optfunc = @(Zdefoc) defoc_func(Zdefoc,phi,imDM(:,:,:,indx),H0,Z4);
    % Evaluate using both signs of the defocus (since it's ambiguous) and
    % choose the smaller
    deferr = zeros(1,2);
    signDefocus = [1 -1];
    for signIndex = 1:2
      deferr(signIndex) = optfunc(signDefocus(signIndex)*Zdefocus(indx));
    end
    [thiserr,signIndex] = min(deferr);
    Zdefoc = signDefocus(signIndex)*Zdefocus(indx);
%     % Just evaluate at a discrete set of values (don't need high precision)
%     n_defoc = length(Zdefoc_values);
%     err_defoc = zeros(1,n_defoc);
%     for defocIndex = 1:n_defoc
%       err_defoc(defocIndex) = optfunc(Zdefoc_values(defocIndex));
%     end
%     [thiserr,minIndex] = min(err_defoc);
%     Zdefoc = Zdefoc_values(minIndex);
    %[Zdefoc,thiserr,exitflag,output] = fminsearch(optfunc,0);
    if options.debug
      %fprintf('Zdefoc = %g (%d evaluations)\n',Zdefoc,output.funcCount);
      fprintf('Zdefoc = %g\n',Zdefoc);
      shg
      for k = 1:K
        phi(:,:,k) = phi(:,:,k) + Zdefoc * Z4;
      end
      [thiserr,g,obj] = pdpenalty(phi,imDM(:,:,:,indx),H0);
      [Hk,sk,imc] = pd_forward_model_2d(phi,H0,obj);
      for k = 1:K
        subplot(1,2,1);
        imshowsc(imDM(:,:,k,indx));
        title(['Actuator ' num2str(act_num)])
        subplot(1,2,2)
        imshowsc(imc(:,:,k));
        drawnow
      end
      pause
    end
    err = err + thiserr;
  end
end

function [err,phi] = defoc_func(Zdefoc,phi,imDM,H0,Z4)
  K = size(phi,3);
  for k = 1:K
    phi(:,:,k) = phi(:,:,k) + Zdefoc * Z4;
  end
  err = pdpenalty(phi,imDM,H0);
end
  