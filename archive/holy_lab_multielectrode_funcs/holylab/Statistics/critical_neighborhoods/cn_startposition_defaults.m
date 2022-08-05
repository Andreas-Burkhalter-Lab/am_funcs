function [params,x0,C] = cn_startposition_defaults(params,varargin)
% Syntax:
%   [params,x0,C0] = cn_startposition_defaults(params,x0)
%   [params,x0,C0] = cn_startposition_defaults(params,x0,C0)
%   [params,x0,C0] = cn_startposition_defaults(params,info)
  
  C = [];
  if isstruct(varargin{1})
    % We initialize with x0 and C0 rather than x and C. This is so we can
    % properly support the setting of params.updateCovariance without
    % making the input syntax complicated ("...uses C0 if updateCovariance is
    % false, otherwise uses the field C..."). The apparent negative to this
    % is that the first iteration climbs back up to the peak, which seems
    % redundant. However, the bonus is that it will automatically terminate
    % the next time the peak is visited, saving iterations.
    x0 = varargin{1}.x0;
    C = varargin{1}.C0;
  else
    x0 = varargin{1};
    if length(varargin) > 1
      C = varargin{2};
    end
  end
  % Initialize covariance & covarianceModel, if necessary
  if ~isfield(params,'covarianceModel') || isempty(params.covarianceModel)
    if isempty(C)
      params.covarianceModel = 'isotropic';
    else
      if (numel(C) == 1)
        params.covarianceModel = 'isotropic';
      elseif (numel(C) == d)
        params.covarianceModel = 'diagonal';
      elseif (numel(C) == d*d)
        params.covarianceModel = 'full';
      else
        error('Can''t determine which covariance model is being used');
      end
    end
  end
  if isempty(C)
    switch params.covarianceModel
      case 'isotropic'
        C = 1;
      case 'diagonal'
        C = ones(d,1);
      case 'full'
        C = eye(d,d);
      otherwise
        error('Covariance model not recognized');
    end
  end