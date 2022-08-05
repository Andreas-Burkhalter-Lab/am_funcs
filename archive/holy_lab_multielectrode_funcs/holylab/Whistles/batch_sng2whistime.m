% BATCH_SNG2WHISTIME: a script to calculate .whistime files from .sng
% files
% This loops over all the .sng files in a given directory. No arguments
% are needed to call this function.
path
p = whistimesdefaults
sngfls = dirbyname('*.sng');
for i = 1:length(sngfls)
  [pathstr,basename] = fileparts(sngfls{i});
  newfilename = [pathstr basename '.whistime'];
  if ~exist(newfilename,'file')
    wt = whistimes(sngfls{i},p)';
    save(newfilename,'wt','-ascii','-double');
  end
  % This matrix is loadable with wt=load(filename)
end
