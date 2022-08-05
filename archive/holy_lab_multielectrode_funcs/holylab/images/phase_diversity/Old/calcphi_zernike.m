function fout = calcphi_zernike(s)
% CALCPHI_ZERNIKE: manage phases from Zernike coefficients
% 
% This is a utility function for calcphi2d, returning the conversion
% functions param2phi (converting parameters to phases) and gphi2gparam
% (converting derivatives with respect to phases into derivatives with
% respect to parameters).
%
% Syntax:
%   f = calcphi_zernike(s)
% where s may have the following fields:
%   rho,theta: pupil coordinates
%   Zindex: a cell array of length K, each entry a vector of Zernike
%     coefficient indices for each diversity image
%   Zcoefs: a cell array of the same configuration as Zindex, giving
%     the coefficient values
%   map: a cell array of the same configuration as Zindex, where
%       each element is a vector containing the index of a parameter
%       vector used for optimization.  0 indicates that the parameter is
%       fixed (not optimized). This notation allows you to "tie"
%       particular coefficients together across diversity images (or you
%       can make them all separate, it's up to you).
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
%
% On output, the fields of f are set to contain the two conversion functions
% needed for CALCPHI2D.  Also set is a field, "extract", which pulls out a
% parameter vector from the initial settings in Zcoefs.
% 
% See also: CALCPHI2D, CALCPHI_LINEARSUMS.
  
% Copyright 2009 by Timothy E. Holy

  persistent Zval ZvalIndex

  %if isfield(s,'A')
  %  s_linear = calcphi_linearsums(s.A,s.offsets);
  %end
  %s.linear = s_linear;
  fout.param2phi = @(p) zern_expand(s,p);
  fout.gphi2gparam = @(grad) zern_contract(grad,s);
  fout.extract = zern_extract(s);

  function phik = zern_expand(s,p,rho,theta)
    % note syntax allows one to provide supplementary rho, theta in case
    % the pupil ever changes?
    Zvalid = true;
    if (nargin > 2)
      Zvalid = false;
    else
      rho = s.rho;
      theta = s.theta;
    end
    framesz = size(rho);

    % Stuff the values of p into cell array of Zernike coefficients
    mapc = s.map;
    K = length(mapc);
    Zcc = s.Zcoefs;
    %map = cat(2,mapc{:});
    for indx2 = 1:K
      thismap = mapc{indx2};
      thismapI = thismap(thismap > 0);
      Zcc{indx2}(thismap>0) = reshape(p(thismapI),size(thismapI));
    end
    
    % If necessary, compute values of Zernike polynomials over the pupil
    uZindex = unique(cat(2,s.Zindex{:}));
    if (~Zvalid || isempty(ZvalIndex) || ~isequal(uZindex,ZvalIndex) || ~isequal(size(Zval(:,:,1)),size(rho)))
      % 
      ZvalIndex = uZindex;
      Zval = zernike_values_normalized(rho,theta,ZvalIndex);
      Zval = cast(Zval,class(rho));
    end
    
    % Compute the phases
    phiksz = [framesz K];
    phik = zeros(phiksz,class(rho));
    for indx = 1:K
      Zindex = findainb(s.Zindex{indx},ZvalIndex);
      phik(:,:,indx) = Zcoefs2phi(Zcc{indx},Zval(:,:,Zindex));
    end
  end
  
  function gradZp = zern_contract(grad,s)
    map = cat(2,s.map{:});
    Zindex = cat(2,s.Zindex{:});
    K = length(s.Zindex);
    Kmap = cell(1,K);
    for indx = 1:K
      Kmap{indx} = repmat(indx,[1 length(s.Zindex{indx})]);
    end
    Kmap = cat(2,Kmap{:});
    clabel = agglabel(map+1); % +1 to handle the 0 map
    clabel(1) = [];
    n = max(map);
    gradZp = zeros(n,1);
    for indx = 1:n
      thisindx = clabel{indx};
      for indx2 = 1:length(thisindx)
        thisindx2 = thisindx(indx2);
        % Project onto Zernike polynomials
        tmp = grad(:,:,Kmap(thisindx2)) .* Zval(:,:,ZvalIndex == Zindex(thisindx2));
        gradZp(indx) = gradZp(indx) + sum(tmp(:));
      end
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

function p = zern_extract(s)
  Zc = cat(2,s.Zcoefs{:});
  map = cat(2,s.map{:});
  p = nan(1,max(map)+1);
  p(map+1) = Zc; % +1 to handle the 0 index
  p(1) = []; % to clear the 0 index
end
