function msprof_nnmf_pairwise_plot(basename,sample1,sample2,options)
  if (nargin < 4)
    options = struct;
  end
  options = default(options,'log',true);
  nnmf_file = [basename '_nnmf.mat'];
  filelist_file = [basename '_filelist.mat'];
  s = load(nnmf_file);
  itmp = load(filelist_file,'isotopeinfo');
  i1 = mnpp_sum(s.nnmfresults,sample1);
  i2 = mnpp_sum(s.nnmfresults,sample2);
  figure;
  hax = gca;
  ops = struct('plotfunc',@(pk) mnpp_plot(pk,sample1,sample2,s.file,s.tc,itmp.isotopeinfo),...
    'hax',hax);
  mdexplore([i1;i2],s.nnmfresults,ops);
  if options.log
    set(hax,'XScale','log','YScale','log');
  end
  
function s = mnpp_sum(nnmfresults,indx)
  s = arrayfun(@(p) sum(p.iProf{indx}),nnmfresults);
  
function mnpp_plot(pk,sample1,sample2,file,tc,isotopeinfo)
  figure
  line(tc{sample1},pk.iProf{sample1},'Color','b');
  line(tc{sample2},pk.iProf{sample2},'Color','r');
  xlabel('Elution time (min)')
  title(sprintf('m/z = %g, charge = %d',pk.mzMeanCorrected,pk.charge))
  nIstrc = cell(1,length(isotopeinfo));
  for k = 1:length(isotopeinfo)
    if (k < length(isotopeinfo))
      fmtstr = '%s: %g, ';
    else
      fmtstr = '%s: %g\n';
    end
    nIstrc{k} = sprintf(fmtstr,isotopeinfo(k).symbol,pk.nI(k));
  end
  nIstr = cat(2,nIstrc{:});
  fprintf(nIstr);
  [~,fn1] = fileparts(file{sample1});
  [~,fn2] = fileparts(file{sample2});
  legend(fn1,fn2);
  formula_gui(pk.mzMeanCorrected*pk.charge,isotopeinfo,pk.nI,pk.nIcov);
  
  