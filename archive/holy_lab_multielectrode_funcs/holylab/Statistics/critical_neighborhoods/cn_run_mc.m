function [T2iso,T2diag,T2full] = run_mcT2(saveflag)
  % Run Monte Carlo simulations of Hotelling's T^2 for many different situations
  % Syntax:
  %   [T2iso,T2diag,T2full] = run_mcT2
  %   [T2iso,T2diag,T2full] = run_mcT2(true)
  % The second version saves the results to .mat files in the current
  % directory.
  %
  % See also: cn_mclinesearch, cn_T2thresh_frommc.
  
  % Copyright 2011-2012 by Timothy E. Holy
  
  %% Initialize
  % If we are saving, let's first test for permission before we go to the
  % effort to do a lot of calculation
  if saveflag
    x = 5;
    save testfile x
    delete testfile
  end

  % Set up the parameter lists. Modify these with care, the computation
  % times can get long!
  dlist = unique(round(logspace(0,2,21)));
  dlist = dlist(dlist <= 20);
  nlist = unique(round(logspace(0,3,31)));
  pvaluelist = [0.05 0.02 0.01 0.005 0.002 0.001];
  %pvaluelist = [0.2 0.1];  % good for debugging (fast)
  
  T2base = struct('dlist',dlist,'nlist',nlist,'pvaluelist',pvaluelist);

  %% Do the computations
  % Isotropic
  params = struct('covarianceModel','isotropic','pvalue',pvaluelist(end));
  T2iso = T2base;
  [T2iso.basepoint,T2iso.linesearch] = rmc_run(dlist,nlist,pvaluelist,params);
  if saveflag
    T2data = T2iso;
    save T2iso T2data
  end
  % Diagonal
  params = struct('covarianceModel','diagonal','pvalue',pvaluelist(end));
  T2diag = T2base;
  [T2diag.basepoint,T2diag.linesearch] = rmc_run(dlist,nlist,pvaluelist,params);
  if saveflag
    T2data = T2diag;
    save T2diag T2data
  end
  % Full
  params = struct('covarianceModel','full','pvalue',pvaluelist(end));
  T2full = T2base;
  [T2full.basepoint,T2full.linesearch] = rmc_run(dlist,nlist,pvaluelist,params);
  if saveflag
    T2data = T2full;
    save T2full T2data
  end
end

function [T2cat,T2lscat] = rmc_run(dlist,nlist,pvaluelist,params)
  T2cat = nan(length(dlist),length(nlist),length(pvaluelist));
  T2lscat = nan(length(dlist),length(nlist),length(pvaluelist));
  nsim = ceil(1/min(pvaluelist)^2);
  for i = 1:length(dlist)
%     savestr = struct('d',dlist(i),'nlist',nlist,'covarianceModel',params.covarianceModel,'nsim',nsim);
%     save('cn_mc_testdata','-struct','savestr');
    %[T2,T2ls] = cn_mclinesearch(dlist(i),nlist,params.covarianceModel,nsim,[],0);
    [T2,T2ls] = cn_mclinesearch(dlist(i),nlist,params.covarianceModel,nsim);
    T2s = sort(T2','descend');
    T2lss = sort(T2ls','descend');
    T2thresh = T2s(round(pvaluelist*size(T2s,1)),:);
    T2lsthresh = T2lss(round(pvaluelist*size(T2s,1)),:);
    T2cat(i,:,:) = T2thresh';
    T2lscat(i,:,:) = T2lsthresh';
  end
  T2cat(isnan(T2cat)) = inf;
  T2lscat(isnan(T2lscat)) = inf;
end
