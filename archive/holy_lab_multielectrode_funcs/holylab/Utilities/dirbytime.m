function files = dirbytime(s)
% dirbytime: return filenames matching pattern, sorted by date/time
% filenames = dirbytime(s)
% s is a string giving the name template; pathnames and wildcards work
% correctly (see help for dir).
d = dir(s);
if (isempty(d))
  files = {};
  return;
end
for i = 1:length(d)
  t(i) = datenum(d(i).date);
end
[ts,i] = sort(t);
[files{1:length(d)}] = deal(d.name);
files = files(i);
