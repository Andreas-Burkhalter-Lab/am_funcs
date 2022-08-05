function comp = connected_components(A)
% CONNECTED_COMPONENTS: find the connected components of a graph
% Syntax:
%   comp = connected_components(A)
% where
%   A is a matrix with A(i,j)=1 whenever nodes i & j are connected;
%   comp is a cell array, where each element is the list of indices in a
%     connected component.
  
% Copyright 2005-2010 by Timothy E. Holy
  
  A=A|A'; % Ensure the matrix is symmetric
  A = A + speye(size(A));  % ensure that the diagonal is marked as "connected"
  n_points = size(A,1);
  map = 1:n_points;
  mapOld = [];
  while ~isequal(map,mapOld)
    mapOld = map;
    for i = 1:n_points
      indx = find(A(:,i));
      map(indx) = min(map(indx));
    end
  end
  [~,~,lbl] = unique(map);
  comp = agglabel(lbl);
  return


  AnnzOld = 0;
  
  A=A|A'; % NOTE: to make the matrix symmetric
  
  A = double(A);
  [i,j] = find(A);
  Annz = length(i);
  while (Annz > AnnzOld)
    AnnzOld = Annz;
    for k = 1:Annz;
      A(i(k),:) = A(i(k),:) + A(j(k),:);
    end
    [i,j] = find(A);
    Annz = length(i);
  end
  A = (A > 0);
  uA = unique(A,'rows');
  ncomp = size(uA,1);
  comp = cell(1,ncomp);
  for i = 1:ncomp
    comp{i} = find(uA(i,:));
  end
