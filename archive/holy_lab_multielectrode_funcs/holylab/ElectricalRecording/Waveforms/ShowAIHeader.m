function ShowAIHeader(header)
% ShowAIHeader(header)
% Prints user-relevant header info on screen
% (all the info can of course be accessed directly)
fprintf('%s\t\t%s\n',header.date,header.time)
disp(header.usrhdr)
if (isfield(header,'version') & header.version == 1)
        disp(['Channels: ',header.chstr])
        fprintf('Acquisition time: %f s\n',header.acqtime)
else
        fprintf('Channels:');
    if (isfield(header,'channels'))
        fprintf(' %d',header.channels);
    else
        fprintf(' %d',header.chan);
    end
        fprintf('\n');
        fprintf('Acquisition time: %f s\n',header.nscans/header.scanrate)
end
fprintf('Scan rate: %f Hz\n',header.scanrate)
