function options=progress_bar(options)
% progress_bar: create/update/destroy a progress bar
% syntax:
%    options=progress_bar(options)
% pre:
%    options: a struct whose valid fields are:
%        progress: if ==-1 or >= max, hide the progress bar
%                  else create/update the progress bar
%        max: the progress' value is [0, max]
%        caption='%d %% done': the caption of the progress bar
%        what="Please wait ...': the description to the progress bar
%        tex_mode=0: use TeX interpreter on "what" if true.
% note:
%    this is a wrapper to progress.m.
% usage:
%    ex. 1:
%        for(idx=1:nFramesToProcess)
%          % do something here
%          progress_bar(struct('progress', idx, 'max', nFramesToProcess));
%        end % for, 
% 
%    ex. 2:
%        for(idx=1:nFramesToProcess)
%          % do sth here
%          if(mod(idx,10)==1 || idx==nFramesToProcess)
%             progress_bar(struct('progress', idx, 'max', nFramesToProcess, 'what', ['converting ' filename '...']));
%          end % if, need to show progress every 10 frames or the last frame
%        end % for, 
% see:
%    progress.m
% 
   if(~isfield(options, 'progress') || ~isfield(options, 'max') )
      error('progress and max must be provided in the argument options');
   end
   
   if(~isfield(options, 'handle'))
      f = findobj(allchild(0),'flat','Tag','TMWWaitbar');
      if(~isempty(f))
	 options.handle=f(1);
      end
   end
   
   if(options.progress>=options.max) 
      options.progress=-1; 
   end % if, 100%, then close progress bar automatically
   
   options=progress(options);
   
