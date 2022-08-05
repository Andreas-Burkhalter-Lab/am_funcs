function files = dirbyname(s)
% DIRBYNAME: return filenames matching pattern, sorted alphabetically
% filenames = dirbyname(s)
% s is a string giving the name template; pathnames and wildcards work
% correctly (see help for dir).
d = dir(s);
if (isempty(d))
  files = {};
  return;
end
[files{1:length(d)}] = deal(d.name);
