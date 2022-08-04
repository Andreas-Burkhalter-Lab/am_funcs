function mouseAccelerationCheck
%mouseAccelerationCheck: pause code execution until mouse acceleration is
%turned off (for measuring mouse movements).
%%% These parameters produce linear mouse movement:cursor movement ratio
%%% only after using the MarkC_Windows MouseFix and as long as the cursor
%%% is prevented from moving to any screen edge.
%%% last updated 9/17/15 on stim comp

MouseThreshold1 = winqueryreg('HKEY_CURRENT_USER', 'Control Panel\mouse','MouseThreshold1');
MouseThreshold2 = winqueryreg('HKEY_CURRENT_USER', 'Control Panel\mouse','MouseThreshold2');
MouseSensitivity = winqueryreg('HKEY_CURRENT_USER', 'Control Panel\mouse','MouseSensitivity');

while ~strcmp(MouseThreshold1,'0') || ~strcmp(MouseThreshold2,'0') % if enhaned pointer precision is still on
    input('Enhanced pointer precision is on. Turn off enhanced pointer precision then press Enter.')
    MouseThreshold1 = winqueryreg('HKEY_CURRENT_USER', 'Control Panel\mouse','MouseThreshold1');
    MouseThreshold2 = winqueryreg('HKEY_CURRENT_USER', 'Control Panel\mouse','MouseThreshold2');
end

while ~strcmp(MouseSensitivity,'10') % if mouse sensitivity is not at the correct value for linear readout
    input(['Mouse sensitivity is set to ' MouseSensitivity...
        '. Set mouse sensitivity to 10 then press Enter.'])
    MouseSensitivity = winqueryreg('HKEY_CURRENT_USER', 'Control Panel\mouse','MouseSensitivity');
end