function intensity_batch
   headerFiles={};
   roidefFiles={};
   
   while(1)
      file_filter='*.roidef';
      [filename, pathname] = uigetfile(file_filter, 'Pick a ROI definition file');
      if(filename==0)
	 return;
      end
      roidefFiles{end+1} =fullfile(pathname, filename);
   
      file_filter='*';
      [filename, pathname] = uigetfile(file_filter, 'Pick an Imagine header file');
      if(filename==0)
	 return;
      end
      headerFiles{end+1} =fullfile(pathname, filename);
      
      msgCont='continue selecting files';
      msgDone='done with selecting files';
      usrResponse=questdlg('Done with this round. What next?', ...
			   'batch calculating stack ROI intensities', ...
			   msgCont, msgDone, msgCont);

      switch(usrResponse)
	 case msgCont
	    continue;
	 case msgDone
	    break;
	 otherwise
	    continue; % if user close questdlg instead of click one of the
		      % two buttons
      end % switch
   end % while,
   
   for idx=1:length(roidefFiles)
      disp(['now processing file #' num2str(idx) ' of ' num2str(length(roidefFiles))]);
      
      inten=calc_stack_roi_inten(headerFiles{idx}, roidefFiles{idx});
      intenFile=replace_extension(roidefFiles{idx}, '.intensity');
      save(intenFile, 'inten', '-mat');
   end
