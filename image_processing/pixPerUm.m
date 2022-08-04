%%% get pixels per micron from zoom value
% updated 2020/11/17
% 'zoom' argin specifies the zoom of the lens used to acquire to picture....
%       if from ab lab epifl. scope starting around May 2018, must be one of these values: 2, 4, 10, 20, 40, 100; 
%           -assuming 1392x1040 pix image for ab lab scope
%       if from ab lab epifl. scope BEFORE MAY 2018, must be one of these values: 1.6, 3.2, 8, 16, 32; 
%           -assuming 1392x1040 pix image for ab lab scope
%       if from wucci olympus confocal scope, must be one of: 4, 10, 20;
%           -assuming 1024x1024 pix image for wucci confocal scope
%       this value will be used to convert pixels into distances in microns
% 'scope'  -if not from ablab epi scope on counter, arg 2 = 'macro','conf' (confocal), 'macro1280x1024'  

function pixelsPerUM = pixPerUm(zoom,scope)

if ~exist('scope','var') || isempty(scope) || strcmp(scope,'epimicro') %%% AB lab epi scope starting ~may 2018
    switch zoom    %  assuming 1392*1040 pixels
        case 2
            pixelsPerUM = 0.3123;
        case 4
            pixelsPerUM = 0.6208;
        case 10 
            pixelsPerUM = 1.5335;
        case 20
            pixelsPerUM = 3.0670;
        case 40
            pixelsPerUM = 6.1339;
        case 100
            pixelsPerUM = 15.3346;
        otherwise
            error('''zoom'' arg must be one of the following: 2, 4, 10, 20, 40, 100 (for scope == ''epimicro'')')
    end
else 
    switch scope
        case 'epimicro_old' %  1392*1040 pixels
            switch zoom
                case 1.6  % ab lab epifl. scope BEFORE MAY 2018
                    pixelsPerUM = 0.2561; 
                case 3.2  % ab lab epifl. scope BEFORE MAY 2018
                    pixelsPerUM = 0.51;
                case 8  % ab lab epifl. scope BEFORE MAY 2018
                    pixelsPerUM = 1.2741;
                case 16  % ab lab epifl. scope BEFORE MAY 2018
                    pixelsPerUM = 2.5501; 
                case 32  % ab lab epifl. scope BEFORE MAY 2018
                    pixelsPerUM = 5.089;
                otherwise	
                    error('''zoom'' arg must be one of the following: 1.6, 3.2, 8, 16, 32 (for scope == ''epimicro_old'')')
            end
        case 'kepecs_thunder' %  2048 x 2048 pixels
            switch zoom
                case 40 % image dimensions = 330.79um x 330.79um (see image metadata)
                    pixelsPerUM = 6.1912;
                otherwise
                     error('''zoom'' arg must be one of the following: 40 (for scope == ''kepecs_thunder'')')
            end
        case 'macro1280x1024' % C:\Users\AM\Documents\sections\from_ABlab\from_AB\09060m2AChR\Left
            switch zoom
                case 2
                    pixelsPerUM = 0.3006;
                case 4 
                    pixelsPerUM = 0.6012; % inferred
                otherwise
                    error('')
            end
        case 'jietalfig3b' % ji et al. fig3b scale bar superimposed on 09060m2AChR 1280x1024 image labeled as x4 zoom
            switch zoom
                case 4
                    pixelsPerUM = 0.5341; % scale bar too big? - underestimates patch density
            end
        case 'jietalfig3b_fullv1' % ji et al. fig3b scale bar superimposed on 09060m2AChR 1792x2125 image showing full v1 labeled as x4 zoom
            switch zoom
                case 4
                    pixelsPerUM = 0.6676; % scale bar too big? - underestimates patch density
            end
        case 'macro' % ablab surgical macro scope, 1392*1040
            switch zoom
                case 1.25
                    pixelsPerUM =  0.1124; % from C:\Users\AM\Documents\sections\Scale 100um_10um_Macro_Scope x1_25
                case 2.5
                    pixelsPerUM =  0.2282;   % from C:\Users\AM\Documents\sections\Scale 100um_10um_Macro_Scope x2_5
            otherwise
                    error('''zoom'' arg must be one of: 1.25, 1.8, 2.5, 3.2 (macro)')
            end
        case 'confocal'
            switch zoom %%% conversions for ab lab scope (measured empirically)
                case 4 % wucci olympus confocal scope
                    pixelsPerUM = 0.32; 
                case 10 % wucci olympus confocal scope
                    pixelsPerUM = 0.7911;
                case 20 % wucci olympus confocal scope
                    pixelsPerUM = 1.5911; % averaged from x4 and x10, not directly measured  
             otherwise
                    error('''zoom'' arg must be one of the following: 4, 10, 20 (WUCCI confocal).')
            end
    end
end
