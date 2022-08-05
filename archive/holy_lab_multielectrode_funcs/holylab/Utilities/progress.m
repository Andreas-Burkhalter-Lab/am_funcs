function options=progress(options)
%  progress: create/update/destroy a progress bar
%  syntax:
%     options=progress(options)
%  pre:
%     options=[]: a struct whose valid fields are:
%        handle: if no such field create a new progress bar;
%                else if is a handle, update the related progress bar;
%                else do nothing
%        progress=0: if ==-1, hide the progress bar
%                    else update the progress bar
%        max=100: the progress' value is [0, max]
%        caption='%d %% done': the caption of the progress bar
%        what="Please wait ...': the description to the progress bar
%        tex_mode=0: use TeX interpreter on "what" if true.
%  post: 
%     if a field is omitted in the arg, it may be filled in the return value
%     especially the handle.
%  usage:
%     ex. 1:
%        tt=progress; % init
%        tt.progress=20; tt=progress(tt); % 20% done
%        tt.progress=-1; tt=progress(tt); % fini
%     ex. 2: 
%        tt=progress(struct('what', ['converting file ' filename ', \nplease wait ...'], 'max', 1024));
%        tt.progress=256; tt=progress(tt); % update progress bar
%        tt.progress=-1;  tt=progress(tt); % finished
% 
   if(nargin==0) options=[]; end
   
   if(~isfield(options, 'handle') || options.handle==-1 )
      options.handle=waitbar(0,'', 'name', '');
   end % if, no handle provided or invalid handle
   
   if(~ishandle(options.handle)) return; end
   
   if(~isfield(options, 'progress'))
      options.progress=0;
   end

   if(options.progress==-1)
      close(options.handle);
      options.handle=-1;
      return ;
   end % if, need to close the progress bar
   
   if(~isfield(options, 'max'))
      options.max=100; % the minimum value is always 0
   end

   if(~isfield(options, 'caption'))
      options.caption=[' %d %% done'];
   end
   
   if(~isfield(options, 'what'))
      options.what='Please wait ...';
   end
   
   if(~isfield(options, 'tex_mode'))
      options.tex_mode=0;
   end

   ttRatio=options.progress/options.max;
   set(options.handle, 'name',   sprintf(options.caption, round(ttRatio*100)) );
   waitbar(ttRatio, options.handle, sprintf(options.what, round(ttRatio*100)) );
   htitle = get(get(options.handle,'Children'),'Title');
   if(options.tex_mode)
      set(htitle,'Interpreter','tex');
   else
      set(htitle,'Interpreter','none');
   end
   