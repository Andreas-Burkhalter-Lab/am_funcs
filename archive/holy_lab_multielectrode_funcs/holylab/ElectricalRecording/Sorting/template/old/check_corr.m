   % fittingResultFile: a *.fit file
   fittingResultFile='cycle2_berry.fit';
   
   fittingResultVars=load(fittingResultFile, '-mat');
   fileToFit=fittingResultVars.fileToFit;
   templateFile=fittingResultVars.templateFile;
   fitting=fittingResultVars.fitting;
   spikeCollection=fittingResultVars.spikeCollection;
   T=fittingResultVars.T;
   ShiftMatrix=fittingResultVars.ShiftMatrix;
   
   templateVars=load(templateFile, '-mat');
   templates=templateVars.templates;
   channels=templateVars.channels;
   thresh=templateVars.thresh;
   medv=templateVars.medv;

  
   nTemplates=length(fitting);
   
   
   curTemplateIdx=35; % template that user is interested.
   
   templatesMat=cell2mat(templates);
   curTemplate=templates{curTemplateIdx};
   curTemplate=repmat(curTemplate, 1, nTemplates);
   distance=sum((templatesMat-curTemplate).^2).^0.5;
   [tt, closestTemplateIndices]=sort(distance);
   closestTemplateIndices=closestTemplateIndices(1:4); % only pick closest three. 
   % NOTE: the first is cur template itself
   spikeTimes=cell(1, length(closestTemplateIndices));
   for idx=1:length(spikeTimes)
      spikeTimes{idx}=[fitting{closestTemplateIndices(idx)}.shiftedTime];
   end % for, each template to compare

   titles=num2cell(closestTemplateIndices);
   [titles{:}]=foreach_g(@num2str, titles{:});
   plot_all_corr(spikeTimes, struct('titles', {titles}));

      