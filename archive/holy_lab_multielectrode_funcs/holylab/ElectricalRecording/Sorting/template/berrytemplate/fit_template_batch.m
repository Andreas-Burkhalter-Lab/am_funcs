function fit_template_batch(templateMerecFile, filesToFit, optionsForFit)
% fit_template_batch:
% PRE: 
%    templateMerecFile: the .merec that is used to generate the templates.
%    filesToFit: a string or cell array of strings
%    optionsForFit: the options passed to fit_template()

   if(nargin==2)
      optionsForFit=struct;
   end
   
   if(~iscell(filesToFit))
      filesToFit={filesToFit};
   end
   
   for fileIdx=1:length(filesToFit)
      curFileToFit=filesToFit{fileIdx};
      for part=1:2
         fineClusterFile=replace_extension(templateMerecFile, ...
            sprintf('_part%d.fine_cluster', part));
         fileToSave=replace_extension(curFileToFit, ...
            sprintf('_part%d.fit', part));
         optionsForFit.fileToSave=fileToSave;
         fit_template(fineClusterFile, curFileToFit, optionsForFit);
      end
   end
