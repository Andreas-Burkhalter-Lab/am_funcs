%%% umPerPix (inverse of pixPerUm)
%%%%% umPerPixel = umPerPix(zoom,scope)
% convert microns to pixels in length
% 'zoom' argin specifies the zoom of the lens used to acquire to picture....
%       if from ab lab epifl. scope, must be one of these values: 1.6, 3.2, 8, 16, 32; 
%           -assuming 1392x1040 pix image for ab lab scope
%       if from wucci olympus confocal scope, must be one of: 4, 10, 20;
%           -assuming 1024x1024 pix image for wucci confocal scope
%       this value will be used to convert pixels into distances in microns
% 'scope'  -if not from ablab epi scope on counter (after ~may 2018), arg 2 = 'epimicro_old','macro','conf' (confocal), 'macro1280x1024'  

function umPerPixel = umPerPix(zoom,scope)

if ~exist('scope','var') || isempty(scope) %  1392*1040
    scope = [];
end

umPerPixel = 1 ./ pixPerUm(zoom,scope);