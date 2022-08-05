function isotopeinfo = isotopes(options)
% ISOTOPES: load NIST information about masses and abundance of isotopes
% This is a wrapper for isoDalton_get_isotope_info by Ross Snider.
% Syntax:
%   isotopeinfo = isotopes;
%   isotopeinfo = isotopes(options);
% where
%   options can contain the following fields:
%     thresh (default 1e-5): isotopes with lower abundance than this are
%       discarded
%     element_list (default {'H','C','O','S','N','P','Na'}): the list of
%       elements for which isotope information will be retrieved
%     n_min (default [0 0 0 0 0 0 0]): for use in mass2formula, the minimum
%       number of each of the elements in any formula
%     n_max (default [Inf Inf Inf 2 2 2 1]): for use in mass2formula, the
%       maximum number of the corresponding element in any formula.
% and
%   isotopeinfo is a structure array, where each element contains
%     information about a single element (see output for details).
%
% See also: MASS2FORMULA.

% Copyright 2009 by Timothy E. Holy

  if (nargin < 1)
    options = struct;
  end
  options = default(options,'thresh',1e-5,...
    'element_list',{'H','C','O','S','N','P','Na'},...
    'n_min',[0 0 0 0 0 0 0],...
    'n_max',[Inf Inf Inf 2 2 2 1]);

  elements = isoDalton_get_isotope_info;
  symbol = cellfun(@(p) p.symbol,elements,'UniformOutput',false);
  indx = findainb(options.element_list,symbol);
  for i = 1:length(indx)
    thiselement = elements{indx(i)}.isotopes;
    abundance = cellfun(@(p) p.isotopic_composition,thiselement);
    mass = cellfun(@(p) p.relative_atomic_mass,thiselement);
    [abundance,sortOrder] = sort(abundance,'descend');
    mass = mass(sortOrder);
    keepflag = abundance > options.thresh*abundance(1);
    abundance = abundance(keepflag);
    mass = mass(keepflag);
    isotopeinfo(i) = struct('symbol',options.element_list{i},'base_mass',mass(1),'isotope_masses',mass,'isotope_abundance',abundance,'n_max',options.n_max(i),'n_min',options.n_min(i));
  end
end
