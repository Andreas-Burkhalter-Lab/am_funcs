function loadindex = ephystofetch(ephysin,fieldname)
% EPHYSTOFETCH: check fields and return the ones needing retrieval
% loadindex = ephystofetch(ephysin,fieldname)
% where
%   ephysin is a structure array of type EPHYS.
%   fieldname is the name of the field to examine
% and
%   loadindex is an index vector indicating the components
%     ephysin(loadindex) which require data to be fetched from disk.
%
% See also: EPHYS, EPHYSFETCH.
nstruct = prod(size(ephysin));
if ~isfield(ephysin,fieldname)
  loadindex = 1:nstruct;
  return;
end
loadii = zeros(1,nstruct);
for i = 1:nstruct
  loadii(i) = isempty(ephysin(i).(fieldname) );
end
loadindex = find(loadii);
