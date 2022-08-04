function code = photoshoot
% photoshoot   Code for the ViRMEn experiment photoshoot.
%   code = photoshoot   Returns handles to the functions that ViRMEn
%   executes during engine initialization, runtime and termination.


% Begin header code - DO NOT EDIT
code.initialization = @initializationCodeFun;
code.runtime = @runtimeCodeFun;
code.termination = @terminationCodeFun;
% End header code - DO NOT EDIT



% --- INITIALIZATION code: executes before the ViRMEN engine starts.
function vr = initializationCodeFun(vr)

% vr.window.WindowState = System.Windows.Forms.FormWindowState.Normal;
% vr.window.Width = 200;
% vr.window.Height = 500;
% vr.oglControl.Width = vr.window.Width;
% vr.oglControl.Height = vr.window.Height;
% aspectRatio = double(vr.oglControl.Size.Width)/double(vr.oglControl.Size.Height);
% cf.oglControl.SetBounds(-aspectRatio,-1,aspectRatio,1);
% GL.Ortho(-aspectRatio,aspectRatio,-1,1,-1e3,0); % orthographic projection

% --- RUNTIME code: executes on every iteration of the ViRMEn engine.
function vr = runtimeCodeFun(vr)



% --- TERMINATION code: executes after the ViRMEn engine stops.
function vr = terminationCodeFun(vr)
