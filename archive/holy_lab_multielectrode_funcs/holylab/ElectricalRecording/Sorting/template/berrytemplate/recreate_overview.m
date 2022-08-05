function recreate_overview(fittingResultFiles, overviewFilename)
% re-generate overview file
% SYNTAX:
%     recreate_overview(fittingResultFiles, overviewFilename)
% PRE:
%    fittingResultFiles: a cell array of .fit filenames
%    overviewFilename: the overview filename to save.

   if(nargin~=2)
      error('Please specify two inputs');
   end
   
   ssnpFiles=cell(1, length(fittingResultFiles));
   for idx=1:length(ssnpFiles)
      ssnpFiles{idx}=replace_extension(fittingResultFiles{idx}, '_fake.ssnp');
   end
   
   overview.sorthead = snipfile2sortheader(ssnpFiles);
   overview.options = []; % TODO: as a place holder
   overview.isFake  = 1;
   
   save(overviewFilename, '-struct', 'overview', '-mat');
