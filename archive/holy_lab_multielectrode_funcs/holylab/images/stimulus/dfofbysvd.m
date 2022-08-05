function R = dfofbysvd(I,trange)
% DFOFBYSVD: compute scalar response using SVD
% Syntax:
%   R = dfofbysvd(I,trange)
% where
%   I is a cell array, one entry per stimulus. Each entry is an array of
%     size n_cells-by-T-by-n_trials, where n_trials is
%     the # of trials and T is the duration of the timerange
%     (diff(trange)+1), containing the response on each trial. Each is
%     expressed as a deltaF/F, where the baseline period ranges from the
%     beginning of the trial's timerange to just before valve onset.
%   trange is a 2-vector, e.g., [-5 10], containing the time interval
%     surrounding valve onset.
% and
%   R is a cell array with n_stim elements, each n_cells-by-n_trials,
%     containing a scalar-valued response for each cell on every trial.
%
% See also: GET_TRIALS_DFOF_CORRECTED.

% Copyright 2010 by Timothy E. Holy

  n_stim = length(I);
  R = cell(1,n_stim);
  % Quantify by SVD on the temporal profiles
  order = [2 1 3];
  tindx = trange(1):trange(2);
  for i = 1:n_stim
    A = I{i};
    % Make successive trials seem like additional cells
    A = permute(A,order);
    Asnip = A(tindx>=0,:,:);
    Asnip(isnan(Asnip)) = 0;
    [u,~,~] = svds(Asnip(:,:),1);
    u = u * sign(mean(u));
    u = max(u,0);
    u = u/sum(u);
    R{i} = squeeze(sum(bsxfun(@times,u,Asnip),1));
  end
