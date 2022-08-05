function optics_out = tune_dm_bruteforce(optics_in,rays_in,targets,minopts)
% TUNE_DM: optimize the parameters of a deformable mirror
% Syntax:
%   optics_out = tune_dm(optics_in,rays,targets,minopts)
% optics_in is a cell array of 2-cells, one of each optical element; the
% first element of every 2-cell is the function handle to the
% "controlling" function and the second the structure that contains the
% parameters describing the element.
% rays is a cell array of ray structure arrays; all the rays in each cell
% are assumed to originate from the same point, so the error penalty is to
% calculate the degree of non-convergence of these rays.
% targets is a n_groups-by-2 matrix, with each row giving the desired
% target point for each group of rays. If empty, it uses the average final
% position of each groups as the target.
% This now tunes the DMs in reverse order: start with the last DM, and then
% re-tune the last two, and then so on.

  optics_out = optics_in;
  if ~iscell(rays_in)
    rays = {rays_in};
  end
  n_rays = length(rays_in);
  % Find the DMs in the optics structure
  n_optics = length(optics_out);
  isDM = false(1,n_optics);
  for optIndex = 1:n_optics
    isDM(optIndex) = isfield(optics_out{optIndex}{2},'DMparams');
  end
  dmIndex_total = find(isDM);
  %for iterIndex = length(dmIndex_total):-1:1
  %for iterIndex = 1:length(dmIndex_total)
    iterIndex = 1;
    dmIndex = dmIndex_total(iterIndex:end);
    %dmIndex = dmIndex_total(1:iterIndex);
    rays = rays_in;
    % Trace the rays up to the first DM, and then discard the previous optics
    % (that way we only need to trace from the point where things
    % start changing)
    if (dmIndex(1) > 1)
      for i = 1:length(rays)
        rays{i} = raytrace(rays{i},optics_out(1:dmIndex(1)-1),false);
        validFlag = [rays{i}.valid];
        rays{i} = rays{i}(validFlag); % keep only valid rays
      end
      optics_out_junked = optics_out(1:dmIndex(1)-1);
      optics_out = optics_out(dmIndex(1):end);
      dmIndex = dmIndex - (dmIndex(1) - 1);
    else
      optics_out_junked = [];
    end
    % Copy over the current DM parameters
    n_dms = length(dmIndex);
    pdmc = cell(1,n_dms); % parameters of DMs as cell array
    lpdm = zeros(1,n_dms); % number of parameters in each DM
    for i = 1:n_dms
      pdmc{i} = optics_out{dmIndex(i)}{2}.DMparams;
      pdm_maxc{i} = optics_out{dmIndex(i)}{2}.DMparams_max;
      lpdm(i) = length(pdmc{i});
    end
    sepIndex = [0 cumsum(lpdm)];
    pdm = cat(2,pdmc{:});
    %tdm_err(pdm)
    pdm_max = cat(2,pdm_maxc{:});
    %pdm = fminsearch(@tdm_err,pdm,minopts);
    %pdm = fminunc(@tdm_err,pdm,minopts);
    pdm = fmincon(@tdm_err,pdm,[],[],[],[],-pdm_max,pdm_max,[],minopts);
    optics_out = [optics_out_junked optics_out];
  %end
  
  function err = tdm_err(p)
    tdm_p2o(p)
    err = 0;
    for i = 1:n_rays
      rf = raytrace(rays{i},optics_out,false);
      if ~all([rf.valid])
        err = inf;
      else
        x0 = [rf.x0];
        if isempty(targets)
          dx = x0 - repmat(mean(x0,2),1,size(x0,2));
        else
          dx = x0 - repmat(targets(:,i),1,length(rf));
        end
        err = err + sum(sum(abs(dx)));
      end
    end
  end

  function tdm_p2o(p)
    % Convert the vector of parameters to the individual DM parameters
    for i = 1:n_dms
      optics_out{dmIndex(i)}{2}.DMparams = ...
        p(sepIndex(i)+1:sepIndex(i+1)); % copy into optics structure
    end
  end
end