function s2 = copyfields(s1,field,s2)
% COPYFIELDS: copy field values from one structure to another
% Syntax:
%   s2 = copyfields(s1,field,s2)
% copies fields listed in the cell array "field" from s1 (if present) and
% adds it to s2.
  
  for i = 1:length(field)
    if isfield(s1,field{i})
      s2.(field{i}) = s1.(field{i});
    end
  end
  