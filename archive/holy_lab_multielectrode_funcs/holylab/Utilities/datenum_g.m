function days=datenum_g(str)
% datenum_g: a wrapper to datenum to parse str in more format
% USAGE: days=datenum_g(str)

   try
      days = datenum(str);
   catch
      days = datenum(str(5:end), 'mmm dd HH:MM:SS yyyy');
   end
   