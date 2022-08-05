function [mass,abundance,formulastring] = isotopologues(formula,isotopeinfo,options)
% ISOTOPOLOGUES: generate the complete list of isotopologues for given molecular formula
%
% Syntax:
%   [mass,abundance,formulastring] = isotopologues(formula,isotopeinfo)
%   [mass,abundance,formulastring] = isotopologues(formula,isotopeinfo,options)
% where
%   formula is a vector containing the number of each element in the
%     formula;
%   isotopeinfo is a structure of the same length as formula;
%   options may have the following fields:
%     thresh (default 0.001): the minimum relative abundance of any
%       isotopologue retained
% and
%   mass is a vector of isotopologue masses
%   abundances is a vector of their corresponding abundances
%   formulastring is a cell array, each listing the isotopic composition.
%
% See also: ISOTOPES.

% Copyright 2009 by Timothy E. Holy

  if (nargin < 3)
    options = struct;
  end
  options = default(options,'thresh',0.001);
  
  % Generate the distribution in each element in the formula
  keepFlag = formula > 0;
  formula = formula(keepFlag);
  isotopeinfo = isotopeinfo(keepFlag);
  n_elements = length(isotopeinfo);
  m = cell(1,n_elements);
  a = cell(1,n_elements);
  f = cell(1,n_elements);
  indx = cell(1,n_elements);
  for i = 1:n_elements
    n_isotopes = length(isotopeinfo(i).isotope_masses);
    if (n_isotopes == 1)
      a{i} = 1;
      m{i} = isotopeinfo(i).base_mass * formula(i);
      f{i} = formula(i);
      indx{i} = 1;
    else
      % Create all combinations of all isotopes having a sum equal to the
      % number of atoms of this element in the formula
      x = cell(1,n_isotopes-1);
      x{1} = 0:formula(i);
      for j = 2:length(x)
        x{j} = x{1};
      end
      if (n_isotopes == 2)
        X = x;
      else
        X = cell(1,n_isotopes-1);
        [X{:}] = ndgrid(x{:});
      end
      X{end+1} = repmat(formula(i),size(X{1}));
      for j = 1:n_isotopes-1
        X{end} = X{end} - X{j};
      end
      Xr = zeros(numel(X{1}),n_isotopes);
      for j = 1:n_isotopes
        Xr(:,j) = X{j}(:);
      end
      % Compute probability of each combination
      pp = isotopeinfo(i).isotope_abundance;
      pp = pp/sum(pp); 
      p = mnpdf_fractional(Xr,pp);
      % Keep the ones that are above threshold
      maxp = max(p);
      keepFlag = p > options.thresh*maxp/n_elements;
      a{i} = p(keepFlag);
      Xr = Xr(keepFlag,:);
      thism = sum(Xr .* repmat(isotopeinfo(i).isotope_masses,sum(keepFlag),1),2);
      m{i} = thism;
      f{i} = Xr;
      indx{i} = 1:size(Xr,1);
    end
  end
  % Combine across elements
  A = cell(1,n_elements);
  M = cell(1,n_elements);
  I = cell(1,n_elements);
  [A{:}] = ndgrid(a{:});
  [M{:}] = ndgrid(m{:});
  [I{:}] = ndgrid(indx{:});
  Atot = A{1};
  Mtot = M{1};
  for j = 2:n_elements
    Atot = Atot .* A{j};
    Mtot = Mtot + M{j};
  end
  Amax = max(Atot(:));
  keepFlag = Atot > options.thresh * Amax;
  mass = Mtot(keepFlag);
  abundance = Atot(keepFlag);
  formulastring = cell(size(mass));
  for j = 1:n_elements
    IkF = I{j}(keepFlag);
    for i = 1:numel(IkF)
      thisN = f{j}(IkF(i),:);
      for k = 1:length(thisN)
        if thisN(k)
          formulastring{i} = [formulastring{i} num2str(round(isotopeinfo(j).isotope_masses(k))) isotopeinfo(j).symbol num2str(thisN(k))];
        end
      end
      formulastring{i} = [formulastring{i} ' '];
    end
  end
  [abundance,sortOrder] = sort(abundance,'descend');
  mass = mass(sortOrder);
  formulastring = formulastring(sortOrder);
end
  
    % Generate via a counter, since most possible outcomes will be of too
    % low abundance to be worth generating    
%     counter = zeros(1,n_isotopes+1);
%     counter(1) = formula(i);
%     indx = 1;
%     while (counter(end) == 0)
%       thisa = mnpdf(counter(1:end-1),isotopeinfo(i).isotope_abundance);
%       if (length(a{i}) > 0)
%         if (thisa < options.thresh*a{i}(1))
%           % 
%           counter(indx) = 0;
%   % Initialize the counter to select the base mass for all atoms
%   basecounter = cell(1,n_elements);
%   for i = 1:n_elements
%     basecounter{i} = [formula(i) zeros(1,length(isotopeinfo.isotope_masses))];
%   end
%   counter = basecounter;
%   counter{n_elements+1} = 0;  % a sentinel for completion
%   mass = [];
%   abundance = [];
%   
%   while (counter{end} == 0)
%     % Compute mass & abundance of this isotopologue
%     thismass = 0;
%     thisab = 1;
%     for i = 1:n_elements
%       thismass = thismass + sum(counter{i} .* isotopeinfo(i).isotope_masses);
%       thisab = thisab * prod(counter{i} .* isotopeinfo(i).isotope_abundance);
%     end
%     if (thisabundance > options.thresh)
%       mass(end+1) = thismass;
%       abundance(end+1) = thisabundance;
%     end
%     % "Increment" the counter
%     i = 1;
%     
  