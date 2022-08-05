function result=check_fit_timing(fittingResultFile)
% result=check_fit_timing(fittingResultFile)
% PRE:
%    fittingResultFile: a *.fit file
% POST:
%    result: 1 if the file need re-fitting.

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

   h=readheader(fileToFit);
   nScans=h.nscans;
   
   nTemplates=length(fitting);
      
   result=1;
   for idx=1:nTemplates
      spikeTimes=[fitting{idx}.shiftedTime];
      if(any(nScans-spikeTimes<47)) % i.e. the distance should >=47. Use 64 for more relax checking
         fprintf('%s need re-fitting\n', fileToFit);
         return;
      end
   end % for, each template

   result=0;
   