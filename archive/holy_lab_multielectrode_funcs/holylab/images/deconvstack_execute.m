function deconvstack_execute(instackfile,psffile,outstackfile)
  smm = stackmm(instackfile);
  load('-mat',psffile);
  m = squeeze(smm(:,:,:,1));
  psf = psf{2};
  psf = reshape(psf,[1 1 length(psf)]);
  %progops = progress_bar(struct('max',nIterations,...
  %                              'what','Deconvolving stack...'));
  img = {m};
  save(outstackfile,'img') % Test for write permission errors
  for i = 1:nIterations
    %progops.progress = i;
    %progress_bar(progops);
    fprintf('%d/%d\n',i,nIterations);
    img = deconvlucy(img,psf,1,0.1);
  end
  save(outstackfile,'img')
