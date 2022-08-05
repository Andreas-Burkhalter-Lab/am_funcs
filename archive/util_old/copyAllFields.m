%COPYALLFIELDS Copies all fields from 'structToCopyFrom' to
%'recipient_struct'.
%  copyallfields(recipient_struct,structToCopyFrom) 
%%%% Identically named fields in recipient_struct will be overwritten by fields 
%%%% from structToCopyFrom.
% last edited 2016/1/6
function recipient_struct_out = copyAllFields(recipient_struct_in,structToCopyFrom) 

recipient_struct_out = recipient_struct_in;

structToCopyFrom_fields = fieldnames(structToCopyFrom);
for f = 1:length(structToCopyFrom_fields)
    recipient_struct_out = setfield(recipient_struct_out, structToCopyFrom_fields{f}, getfield(structToCopyFrom,structToCopyFrom_fields{f}));
end