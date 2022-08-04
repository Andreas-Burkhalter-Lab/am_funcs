function chantest

    handle = cedFunction('SonOpenOldFile', 'Z:\Data\CED\Ben\m1c160r5.smr', 1)
    
    for i = 0:31
        code = cedFunction('SonChanKind', handle, i)
    end
    
    cedFunction('SonCloseFile', handle);