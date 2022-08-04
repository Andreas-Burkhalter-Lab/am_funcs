function code = twoArenas
% twoArenas   Code for the ViRMEn experiment twoArenas.
%   code = twoArenas   Returns handles to the functions that ViRMEn
%   executes during engine initialization, runtime and termination.


% Begin header code - DO NOT EDIT
code.initialization = @initializationCodeFun;
code.runtime = @runtimeCodeFun;
code.termination = @terminationCodeFun;
% End header code - DO NOT EDIT



% --- INITIALIZATION code: executes before the ViRMEN engine starts.
function vr = initializationCodeFun(vr)

vr.cylPos{1} = [55*cosd(20) 55*sind(20)];
vr.cylPos{2} = [-55*cosd(20) 55*sind(20)];


% --- RUNTIME code: executes on every iteration of the ViRMEn engine.
function vr = runtimeCodeFun(vr)

if norm(vr.position(1:2)-vr.cylPos{vr.currentWorld})<20
    numseg = eval(vr.exper.variables.numseg);
    rot = fix(rand*numseg)/numseg*2*pi;
    M = [cos(rot) -sin(rot); sin(rot) cos(rot)];
    vr.worlds{vr.currentWorld}.surface.vertices(1:2,:) = M*vr.worlds{vr.currentWorld}.surface.vertices(1:2,:);
    vr.cylPos{vr.currentWorld} = (M*vr.cylPos{vr.currentWorld}')';
    vr.currentWorld = round(rand)+1;
end


% --- TERMINATION code: executes after the ViRMEn engine stops.
function vr = terminationCodeFun(vr)
