function mutinf = mutual_information(label)
% MUTUAL_INFORMATION: compute the mutual information between labellings
% Syntax:
%   mutinf = mutual_information(label)
% where
%   label is a 2-by-npts matrix of integer labels; these labels determine
%     probability densities in that p_i is equal to the fraction of
%     labels equal to i.
% and
%   mutinf is a scalar equal to the mutual information between
%     probability densities derived from the labels.  

% Copyright 2005 by Timothy E. Holy
  
  [nIter,nPoints] = size(label);
  if (nIter ~= 2)
    error('Input labels must be 2-by-n');
  end
  [upairLabel,tmp,upairIndex] = unique(label','rows');
  [upairList,nList] = agglabel(upairIndex);
  [uLabel1,nLabel1] = agglabel(label(1,:));
  [uLabel2,nLabel2] = agglabel(label(2,:));

  pjoint = nList/nPoints;
  p1 = nLabel1(upairLabel(:,1)')/nPoints;
  p2 = nLabel2(upairLabel(:,2)')/nPoints;
  mutinf = sum(pjoint .* log(pjoint./(p1.*p2)));
  