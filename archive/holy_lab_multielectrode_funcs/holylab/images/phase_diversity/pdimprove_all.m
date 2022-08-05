function [pdparams,err,mu_out] = pdimprove_all(pdparams,mu_in)
% PDIMPROVE_ALL: iterative improvement along distinct coordinate axes for phase diversity
% Syntax:
%  [pdparams_out,err,mu_out] = pdimprove_all(pdparams,mu_in)
% where
%   pdparams is a structure as described in ocpi_pdwrapper;
%   mu_in is a vector of mu values for line minimization (try using ones
%     when in doubt)
% and
%   pdparams_out contains the resulting improved parameters
%   err contains the current mismatch error
%   mu_out contains the values of mu used in the line minimization

  switch_lookup = {'f','sigma','Zcoefs'};
  switch_name = {'fgrad','vgrad','phikgrad'};
  output_name = {'fgrad','sigmagrad','Zgrad'};
  datafields = pdparams.datafields;
  mu_out = mu_in;
  for i = 1:length(datafields)
    indx = strmatch(datafields{i},switch_lookup);
    s = ocpi_pdwrapper(pdparams,struct(switch_name{indx},true));
    g = s.(output_name{indx});
%     if iscell(g)
%       g = cat(2,g{:});  % for Zcoefs derivative
%     end
    g = g(:);
    pdparams.datafields = datafields(i);
    p = ocpi_pdwrapper(pdparams,'extract');
    myfun = @(mu) ocpi_pdwrapper(p - mu*g,pdparams);
    [mu_out(i),err] = lindip(myfun,struct('mu',mu_in(i),'startval',s.val));
    pdparams = ocpi_pdwrapper(p - mu_out(i)*g,pdparams,'insert');
  end
  pdparams.datafields = datafields;
end