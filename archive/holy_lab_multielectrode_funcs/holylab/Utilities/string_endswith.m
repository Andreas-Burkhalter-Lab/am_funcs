function result=string_endswith(str, suffix)
   result=false;
   if(length(str)<length(suffix)) return; end
   result=isequal(str(end-length(suffix)+1:end), suffix);
   
