function Ic = msprof_consolidate_nnmf(nnmfresults)
% msprof_consolidate_nnmf: convert structure to matrix format
%
% This function is useful for analyses such as temporal registration,
% where you are content to discard the m/z information and focus on
% elution of individual "mass channels".
%
% Syntax:
%   Ic = msprof_consolidate_nnmf(nnmfresults)
% where
%   nnmfresults is a structure array of the type saved as output of
%     msprof_nnmf
% and
%   Ic is a 1-by-n_files cell array containing matrices where each row
%     corresponds to a different "mass channel" in nnmfresults.
%
% See also: msprof_temporal_register.
  
% Copyright 2010 by Timothy E. Holy
  
  n_files = length(nnmfresults(1).??);
  n_mz = length(nnmfresults);
  n_samples = cellfun(@length,nnmfresults(1).??);
  Ic = cell(1,n_files);
  for k = 1:n_files
    Ic{k} = zeros(n_mz,n_samples(k));
  end
  for k = 1:n_mz
    for j = 1:n_files
      Ic{k}(j,:) = nnmfresults(k).??{j};
    end
  end
  