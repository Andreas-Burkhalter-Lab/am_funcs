handle = cedFunction('SonOpenOldFile', 'm2c256r2sorted.smr', 1);

if (handle < 0)
    return;
end

code = cedFunction('SonChanKind', handle, 6);

% Maxtime
maxtime = cedFunction('SonChanMaxTime', handle, 6)

% Get Event Data
[num, data] = cedFunction('SonGetEventData', handle, 6, 100000, 0, maxtime);

cedFunction('SonCloseFile', handle);