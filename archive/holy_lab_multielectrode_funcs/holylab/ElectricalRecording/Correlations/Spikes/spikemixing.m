function omega = spikemixing(spikes,params)
% SPIKEMIXING: estimate contamination of clusters from timing info
% Uses the refractoriness of neurons to estimate the purity of clusters
% and the degree to which individual cells are split across cluster
% boundaries.  The mathematical foundation for this is described in
% /usr/lab/doc/refract.pdf.
%
% Syntax:
%   omega = spikemixing(spikes,params)
% where
%   spikes is a nintervals-by-ncells cell array of spike times;
%   params is a structure with the following fields:
%     tref is the duration of the refractory period
%     tp is a 2-vector [start end] giving the period used to estimate the
%       "plateau"; 
%     tpac (optional) is a 2-vector [start end] giving the period used to
%       estimate the "plateau" in autocorrelations, if you want to use
%       something different from tp;
%     tpcc (optional) is a 2-vector [start end] giving the period used to
%       estimate the "plateau" in crosscorrelations, if you want to use
%       something different from tp;
%     plot (optional) if true, causes auto- and cross-correlations to be
%       plotted;
% and
%   omega is an ncells-by-ncells matrix, with the diagonal terms containing
%   the purity for single clusters (Omega_ii), and omega(i,j) contains the
%   mixing fraction omega_i/j. 
%
% See also: AUTOCORRSPIKE, CROSSCORRSPIKE.

  if ~isfield(params,'tpac')
    params.tpac = params.tp;
  end
  if ~isfield(params,'tpcc')
    params.tpcc = params.tp;
  end
  if ~isfield(params,'plot')
    params.plot = 0;
  end
  [nintervals,ncells] = size(spikes);
  if params.plot
    clf
    nbins = 50;
  end
  % Calculate the number of spikes within refractory period,
  % and estimate the number that would be expected from the plateau region
  phi = zeros(ncells,ncells);
  for i = 1:ncells
    tmax = max([params.tref,params.tpac]);
    dt = autocorrspikecat(spikes(:,i),tmax);
    if iscell(dt)
      dt = cat(2,dt{:});
    end
    if params.plot
      subplot(ncells,ncells,(i-1)*ncells + i)
      x = linspace(0,tmax,nbins+1);
      nperbin = histc(dt,x);   % the last element contains points on edge, discard
      xm = (x(1:end-1)+x(2:end))/2;
      % Make the plot symmetric, so it matches the crosscorr plots
      bar([-xm(end:-1:1),xm],[nperbin(end-1:-1:1),nperbin(1:end-1)]);
      set(gca,'XLim',[-1 1]*x(end));
    end
    adt = abs(dt);
    % Number within refractory period
    phi(i,i) = 2*length(find(adt < params.tref));
    % Number in plateau areas
    bigphi(i,i) = 2*params.tref*length(find(adt > params.tpac(1) & ...
      adt <= params.tpac(2)))/diff(params.tpac);
    tmax = max([params.tref,params.tpcc]);
    for j = i+1:ncells
      dt = crosscorrspikecat(spikes(:,i),spikes(:,j),tmax);
      if iscell(dt)
        dt = cat(2,dt{:});
      end
      if params.plot
        subplot(ncells,ncells,(i-1)*ncells + j)
        x = linspace(-tmax,tmax,2*nbins+1);
        nperbin = histc(dt,x);   % the last element contains points on edge, discard
        xm = (x(1:end-1)+x(2:end))/2;
        bar(xm,nperbin(1:end-1));
        set(gca,'XLim',[-1 1]*x(end));
      end
      adt = abs(dt);
      % Number within refractory period
      nref = length(find(adt < params.tref));
      phi(i,j) = nref;
      phi(j,i) = nref;
      % Number in plateau areas
      plat = length(find(adt > params.tpcc(1) & ...
        adt <= params.tpcc(2)))/diff(params.tpcc);
      bigphi(i,j) = params.tref*plat;
      bigphi(j,i) = bigphi(i,j);
    end
  end
  for i = 1:ncells
    omega(i,i) = sqrt(1-phi(i,i)/bigphi(i,i));
  end
  for i = 1:ncells
    for j = i+1:ncells
      omega(i,j) = (bigphi(i,j)-phi(i,j))/(bigphi(j,j)-phi(j,j));
      omega(j,i) = (bigphi(j,i)-phi(j,i))/(bigphi(i,i)-phi(i,i));
    end
  end
  