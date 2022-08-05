function [I,b,Ictrl_recon] = get_trials_dfof_corrected(intensities, onset_list, trange, ctrlindex, options)
% GET_TRIALS_DFOF_CORRECTED: snip temporal traces, subtracting a control timecourse
%
% This obtains temporal traces for each stimulus, and subtracts an estimate
% of the valve artifact.  To use this, at least one of your stimuli must be
% a negative control. The control response is modeled as
%      I(cellnumber, t, trialnumber) = A(cellnumber) * B(trialnumber) * C(t)
% where t is an index corresponding to time within a trial. In subtracting
% this artifact, we take into account that the overall amplitude of the
% artifact can vary, for example with longer gaps between delivery of the
% same valve.
% Note: the fitting of b works reliably only if a minority of cells respond
% to each stimulus. If this is not the case, use the 'linear' mode below.
%
% Syntax:
%   I = get_trials_dfof_corrected(intensities, onset_list, trange, ctrlindex)
%   [I,b,Ictrl] = get_trials_dfof_corrected(intensities, onset_list, trange, ctrlindex)
%   ... = get_trials_dfof_corrected(...,options)
% where
%   intensities is a maxtrix or vector of intensity as a function of time
%     (if supplied as a matrix, is should be of size n_cells-by-n_stacks)
%   onset_list is a cell array, each element containing a vector of valve
%     onset times for each stimulus (see find_stimulus_start)
%   trange is a 2-vector, e.g., [-5 10], specifying the range of time
%     indices (inclusive) relative to the times in onset_list.
%   ctrlindex is the entry of onset_list corresponding to a negative
%     control, whose time course will be subtracted from each response
%     (this is a bit subtle, see above)
%   options may have the following fields:
%     mode: 'independent' or 'linear' (default: 'independent'). In
%       'independent' mode, each valve opening can have its own amplitude B.
%       In 'linear' mode, the amplitude is constrained to be of the form B =
%       b0 + b1*t, where t is the number of stacks since the valve was last
%      openend. (For the first trial, the assumption is that the valve was
%      openend on stack 0.)
% On output,
%   I is a cell array, one entry per stimulus. Each entry is an array of
%     size n_cells-by-T-by-n_trials, where n_trials is
%     the # of trials and T is the duration of the timerange
%     (diff(trange)+1), containing the response on each trial. Each is
%     expressed as a deltaF/F, where the baseline period ranges from the
%     beginning of the trial's timerange to just before valve onset.
%   b is a n_stimuli-by-n_trials matrix, containing the amplitude of the
%     control response subtracted from each trial (in units of the mean
%     control response).
%   Ictrl is the response to the control stimulus.
%
% Note: you can check the dependence of the artifact amplitude on the time
% since last valve opening in the following way:
%   ba = cat(1,b{:});
%   bsnip = ba(:,2:end);  % don't look at the first, cuz don't know time gap
%   dt = zeros(n_stim,n_trials-1); for i = 1:n_stim; dt(i,:) = diff(onset_list{i}); end
%   figure; scatter(dt(:),bsnip(:))
% The time axis is in number of stacks relative to the gap between control
% trials.
%
% See also: GET_TRIALS_DFOF, FIND_STIMULUS_START.

% Copyright 2010 by Timothy E. Holy

  if (nargin < 5)
    options = struct;
  end
  options = default(options,'mode','independent');
  
  n_stim = length(onset_list);
  I = cell(1,n_stim);
  for i = 1:n_stim
    I{i} = get_trials_dfof(intensities,onset_list{i},trange);
  end
  
  % Calculate the control response, so we can subtract it. This response is
  % modeled in the following way:
  %    I(cellnumber,t) = A(cellnumber)*Ir(t)
  % which is optimized by SVD of I with a single component
  Ictrl = I{ctrlindex};
  sz = size(Ictrl);
  Ictrl = permute(Ictrl,[2 1 3]);
  Ictrl = Ictrl(:,:)';
  Ictrl(isnan(Ictrl)) = 0;
  [u,s,v] = svds(Ictrl,1);
  u = reshape(u,sz([1 3]));
  u = mean(u,2);
  Ictrl_recon = u*s*v';
  b = cell(1,n_stim);
  switch options.mode
    case 'independent'
      for i = 1:n_stim
        n_trials = size(I{i},3);
        b{i} = zeros(1,n_trials);
        for j = 1:n_trials
          Itmp = I{i}(:,:,j);
          btmp = robustfit(Ictrl_recon(:),Itmp(:));
          b{i}(j) = btmp(2);
        end
      end
    case 'linear'
      n_trials = cellfun(@(x) size(x,3),I);
      N = sum(n_trials);
      Iall = cell(1,n_stim);
      C0 = cell(1,n_stim);
      C1 = cell(1,N);
      indx = 1;
      for i = 1:n_stim
        % In fitting, we need to omit the first trial, because we don't
        % know how long ago the valve was opened
        dt = diff(onset_list{i});
        Itmp = I{i}(:,:,2:end);
        Iall{i} = Itmp(:);
        C0{i} = repmat(Ictrl_recon(:),n_trials(i)-1,1);
        for j = 2:n_trials(i)
          C1{indx} = dt(j-1)*Ictrl_recon(:);
          indx = indx+1;
        end
      end
      btmp = robustfit([cat(1,C0{:}), cat(1,C1{:})],cat(1,Iall{:}));
      % Now create the actual b values
      for i = 1:n_stim
        dt = diff([0; onset_list{i}]);  % first trial: assume last opened at stack 0
        b{i} = btmp(2) + btmp(3)*dt';
      end
    otherwise
      error('mode not recognized');
  end
  % Perform the subtraction
  for i = 1:n_stim
    btmp = reshape(b{i},[1 1 length(b{i})]);
    Isub = bsxfun(@times,Ictrl_recon,btmp);
    I{i} = I{i} - Isub;
  end
  return
  
  % Compute dot product onto Ictrl
  dp = cell(1,n_stim);
  
  % The amplitude is dependent upon the time since the previous delivery of
  % the same stimulus. 
  
  % Subtract the artifact
  if (nargout < 2)
    for i = 1:n_stim
      I{i} = I{i} - Ictrl_recon;
    end
  else
    b = zeros(n_stim,n_trials);
    for i = 1:n_stim
      % The size of the artifact turns out to depend on how much time has elapsed
      % since that valve was last used. So estimate a single overall scale
      % factor
      %     A = bsxfun(@rdivide,I{i},u);
      %     B = bsxfun(@rdivide,sum(bsxfun(@times,A,V),2),sum(V.^2,2));
      B = bsxfun(@rdivide,sum(bsxfun(@times,I{i},Ictrl_recon),2),sum(Ictrl_recon.^2,2));
      btmp = median(B,1);  % this works if a minority of cells respond to each stimulus
      b(i,:) = squeeze(btmp)';
      I{i} = I{i} - bsxfun(@times,btmp,Ictrl_recon);
    end
  end
  