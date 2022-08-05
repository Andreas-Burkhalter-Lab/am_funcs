function H = wdecomp_amplitude_hessian(w_overlap,t)
% wdecomp_amplitude_hessian: waveform decomposition for the amplitudes
% Syntax:
%   H = wdecomp_amplitude_hessian(w_overlap,t)
% where
%   w_overlap is the set of component overlaps as a function of time
%     shift, see COMPONENT_OVERLAPS.
%   t is the set of event times (a vector; if you have events categorized
%     by type, lump all types together).  This vector must be sorted in
%     increasing order.
% and
%   H is the Hessian matrix used for setting the amplitudes of each
%   waveform component.  The amplitudes A solve the equation
%       H A(:) = B(:)
%   where A(i,j) is the amplitude of the ith component for event j, and
%   B(i,j) is the dot product between the ith component and the voltage
%   waveform for event j.  See the exposition in the preprint for details.
%
% See also: COMPONENT_OVERLAPS, WDECOMP_AMPLITUDE_RHS.
  
% Copyright 2007 by Timothy E. Holy
  
  if ~issorted(t)
    error('The event times must be sorted');
  end
  n_components = size(w_overlap,1);
  l = (size(w_overlap,3)+1)/2;
  n_times = length(t);
  w_overlap = permute(w_overlap,[3 1 2]);
  %comp_t_index = repmat((1:n_components)',1,n_times) + ...
  %    repmat((0:n_times-1)*n_components,n_components,1);
  % Create a lookup that maps (component,event#) to a single vector index
  comp_t_index = reshape((1:n_components*n_times),n_components,n_times);
  % Find the times when there is spike overlap
  [tac,acindx] = autocorrspike(t,l);
  n_pairs = length(tac);
  % Set up storage for making the sparse matrix. We'll use sparse(i,j,s),
  % so create temporary storage for these indices and values.  Thinking in
  % block form (where n_components consist of a block), there is a pair of
  % entries for each overlap, plus the diagonal: hence 2*n_pairs+n_times
  % entries. (But each entry is a block)
  % However, since it's symmetric we can just use the transpose to
  % calculate the lower part.
  % Do the off-diagonal blocks first
  pairI = cell(1,n_pairs);
  pairJ = pairI;
  pairS = pairI;
  for pairIndex = 1:n_pairs
    indx1 = comp_t_index(:,acindx(1,pairIndex));
    indx2 = comp_t_index(:,acindx(2,pairIndex));
    indx1 = repmat(indx1',n_components,1);
    indx2 = repmat(indx2,1,n_components);
    shift = -tac(pairIndex);
    shiftIndex1 = shift+l;
    %shiftIndex2 = -shift+l;
    pairI{pairIndex} = indx1(:);
    pairJ{pairIndex} = indx2(:);
    pairS{pairIndex} = w_overlap(shiftIndex1,:)';
  end
  pairI = cat(2,pairI{:});
  pairJ = cat(2,pairJ{:});
  pairS = cat(2,pairS{:});
  % Do the diagonal blocks
  diagI = cell(1,n_times);
  diagJ = diagI;
  diagS = diagI;
  self_overlap = squeeze(w_overlap(l,:,:));
  self_overlap = self_overlap(:);
  for selfIndex = 1:n_times
    indx = comp_t_index(:,selfIndex);
    indx1 = repmat(indx',n_components,1);
    indx2 = repmat(indx,1,n_components);
    diagI{selfIndex} = indx1(:);
    diagJ{selfIndex} = indx2(:);
    diagS{selfIndex} = self_overlap;
  end
  diagI = cat(2,diagI{:});
  diagJ = cat(2,diagJ{:});
  diagS = cat(2,diagS{:});
  H = sparse([pairI(:);pairJ(:);diagI(:)],[pairJ(:);pairI(:);diagJ(:)],[pairS(:);pairS(:);diagS(:)]);
  