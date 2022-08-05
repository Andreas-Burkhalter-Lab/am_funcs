function result=make_abs_path(apath)
% generate absolute path from a path
% USAGE:
%    result=make_abs_path(apath)
   
   if(isempty(apath))
      error('the arg apath cannot be empty');
   end
   
   if(is_abs_path(apath))
      result=apath;
   else
      if(is_drive_rel_path(apath))
	 error('not implemented yet');
      else
	 result=fullfile(pwd, apath);
      end
   end
   
   
function result=is_abs_path(apath)
   result=0;
   isWindows=ispc;
   if(isWindows)
      if(length(apath)>=3 && apath(2)==':')
	 if(apath(3)=='/' || apath(3)=='\')
	    result=1;
	    return;
	 end
      end
   else
      if(apath(1)=='/')
	 result=1;
	 return;
      end
   end

% test if apath is sth like c:dir   
function result=is_drive_rel_path(apath)
   result=0;
   isWindows=ispc;
   if(isWindows)
      if(~is_abs_path(apath))
	 if(length(apath)>=2 && apath(2)==':')
	    result=1;
	 end
      end
   end
   
