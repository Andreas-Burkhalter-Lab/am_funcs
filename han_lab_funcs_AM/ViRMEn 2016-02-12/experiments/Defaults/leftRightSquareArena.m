function code = leftRightSquareArena
% photoshoot   Code for the ViRMEn experiment leftRightSquareArena.
%   code = leftRightSquareArena Returns handles to the functions that ViRMEn
%   executes during engine initialization, runtime and termination.


% Begin header code - DO NOT EDIT
code.initialization = @initializationCodeFun;
code.runtime = @runtimeCodeFun;
code.termination = @terminationCodeFun;
% End header code - DO NOT EDIT



% --- INITIALIZATION code: executes before the ViRMEN engine starts.
function vr = initializationCodeFun(vr)

% Set current state
vr.currentState = 'cueZone';

% Make south wall invisible and remove associated edge
vr.southWallIndx{1} = vr.worlds{1}.objects.triangles(vr.worlds{1}.objects.indices.southWall,:);
vr.southWallIndx{2} = vr.worlds{1}.objects.triangles(vr.worlds{1}.objects.indices.southWallBottomRim,:);
vr.southWallIndx{3} = vr.worlds{1}.objects.triangles(vr.worlds{1}.objects.indices.southWallTopRim,:);
vr.southWallEdge = vr.worlds{1}.objects.edges(vr.worlds{1}.objects.indices.southWallBottomRim,1);
for ndx = 1:3
    vr.worlds{1}.surface.visible(vr.southWallIndx{ndx}(1):vr.southWallIndx{ndx}(2)) = false;
end
vr.worlds{1}.edges.radius(vr.southWallEdge) = 1e6;


% --- RUNTIME code: executes on every iteration of the ViRMEn engine.
function vr = runtimeCodeFun(vr)

% Check if the animal has entered the main room
if strcmp(vr.currentState,'cueZone') && vr.position(2) > 0
    for ndx = 1:3
        vr.worlds{1}.surface.visible(vr.southWallIndx{ndx}(1):vr.southWallIndx{ndx}(2)) = true;
    end
    vr.worlds{1}.edges.radius(vr.southWallEdge) = 2;
end


% --- TERMINATION code: executes after the ViRMEn engine stops.
function vr = terminationCodeFun(vr)
