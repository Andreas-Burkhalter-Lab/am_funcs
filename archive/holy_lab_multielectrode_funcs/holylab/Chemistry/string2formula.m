function [formula,isotopeinfo] = string2formula(fstr,isotopeinfo)
% string2formula: calculate a numeric "formula" given a chemical formula string
% Syntax:
%   formula = string2formula(fstr,isotopeinfo)  % preferred syntax
%   [formula,isotopeinfo] = string2formula(fstr)
% where
%   fstr is the formula string
%   isotopeinfo is isotope information as loaded by isotopes. You should
%     make sure that all the elements in fstr have data loaded for them in
%     isotopeinfo (unless you're using the non-preferred syntax)
% and
%   formula is a vector containing the mulitiplicities of each element
%     listed in isotopeinfo
%
% Example:
%   isotopeinfo = isotopes;
%   formula = string2formula('C21 H32 O7 S',isotopeinfo)
% yields
%   formula = [32 27 7 1 0 0 0]

% Copyright 2009 by Timothy E. Holy

  fsplit = regexp(fstr,'[a-zA-Z]{1,2}[0-9]*','match');
  elementname = regexp(fsplit,'[a-zA-Z]{1,2}','match');
  abundancestr = regexp(fsplit,'[0-9]*','match');
  n_elements = length(elementname);
  if (length(abundancestr) ~= n_elements)
    error('Parsing error');
  end
  abundance = zeros(1,n_elements);
  for i = 1:n_elements
    if length(elementname{i}) ~= 1
      error('Parsing error');
    else
      elementname{i} = elementname{i}{1};
    end
    if length(abundancestr{i}) > 1
      error('Parsing error');
    else
      if isempty(abundancestr{i})
        abundance(i) = 1;
      else
        abundance(i) = str2double(abundancestr{i}{1});
      end
    end
  end
  if (nargin == 1)
    isotopeinfo = isotopes(struct('element_list',{elementname}));
    formula = abundance;
  else
    formula = zeros(1,length(isotopeinfo));
    try
      formula(findainb(elementname,{isotopeinfo.symbol})) = abundance;
    catch
      error('You are missing one or more elements in isotopeinfo, add them before calling this function');
    end
  end
end
