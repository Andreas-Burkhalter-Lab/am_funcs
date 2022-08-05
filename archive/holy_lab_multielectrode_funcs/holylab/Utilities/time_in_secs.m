function t = time_in_secs(s)
% TIME_IN_SECS: parse time strings, return in units of seconds
% Usage
%   t = time_in_secs(s)
% where s is a string containing a duration, e.g. '102.3ms';
% and t is the time in seconds (e.g., 0.1023)
  
  [t,count,errmsg,nextindex] = sscanf(s,'%g');
  if (count ~= 1)
    error('Error parsing time');
  end
  units = strtrim(s(nextindex:end));
  switch units
   case 's'
   case 'ms'
    t = t/1000;
   case 'min'
    t = t*60;
   case 'hr'
    t = t*3600;
   otherwise
    error(['Time units ' units ' not recognized']);
  end
  