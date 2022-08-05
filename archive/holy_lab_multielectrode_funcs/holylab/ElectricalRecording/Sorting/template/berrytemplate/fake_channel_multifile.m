function fake_channel_multifile(fittingResultFiles, sortResultDir, options)
% fake_channel_multifile(fittingResultFiles, sortResultDir)
% PRE:
%    fittingResultFiles: a string or a cell array of file name strings
%                        (.fit files)
%    sortResultDir: the dir to hold the faked autosort result
%    options: currently only support one field:
%       nClosestTemplates=5: for each fake channel, #templates to consider
% 
% SEE: fake_channel().

   if(nargin==1)
      sortResultDir='sort_fake';
   end
   default_options('nClosestTemplates', 5);
   
   if(~iscell(fittingResultFiles))
      fittingResultFiles={fittingResultFiles};
   end
   
   nFittingResultFiles=length(fittingResultFiles);
   
   % make sure all .fit files use same template file
   templateFiles=cell(1, nFittingResultFiles);
   for fileIndex=1:nFittingResultFiles
      tt=load(fittingResultFiles{fileIndex}, 'templateFile', '-mat');
      templateFiles{fileIndex}=tt.templateFile;
   end
   if(nFittingResultFiles>1)
      if(length(unique(templateFiles))>1)
         error('more than 1 template file are used for the .fit files');
      end
   end
   
   ssnpFiles=cell(1, nFittingResultFiles);
   for fileIndex=1:nFittingResultFiles
      ssnpFiles{fileIndex}=fake_channel(fittingResultFiles{fileIndex}, sortResultDir, ...
                                        struct('onlyGenSsnp', fileIndex>1, ...
                                        'nClosestTemplates', options.nClosestTemplates ...
                                        ) );
   end % for, each .fit file
   
   if(nFittingResultFiles==1)
      return; % skip re-generating overview file b/c it was generated above.
   end
   
   
   % now the overview file
   overview.sorthead = snipfile2sortheader(ssnpFiles);
   overview.options = []; % TODO: as a place holder
   overview.isFake  = 1;
   
   overviewFilename=fullfile(sortResultDir, 'overview.mat');
   save(overviewFilename, '-struct', 'overview', '-mat');
   
   
   
   