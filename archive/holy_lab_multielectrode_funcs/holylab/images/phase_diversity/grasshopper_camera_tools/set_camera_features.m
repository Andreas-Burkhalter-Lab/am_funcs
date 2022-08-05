function src = set_camera_features(src)

set(src,'AutoExposureMode','manual');
set(src,'BrightnessAbsolute',0);
set(src,'BrightnessControl','absolute');
set(src,'GainAbsolute',10);
set(src,'GainControl','absolute');
set(src,'GainMode','manual');
set(src,'Sharpness',1000);
set(src,'SharpnessMode','manual');
set(src,'ShutterAbsolute',0.05);
set(src,'ShutterControl','absolute');
set(src,'ShutterMode','manual')

end