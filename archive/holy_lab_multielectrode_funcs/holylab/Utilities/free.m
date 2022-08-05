function free(handles)
% free: safely delete handles. Like what in C, won't complain if delete a non-handle.
% USAGE: 
%    free(handles)
   delete(handles(ishandle(handles)));
   
