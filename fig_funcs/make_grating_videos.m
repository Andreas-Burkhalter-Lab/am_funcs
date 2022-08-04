%%%% make videos to demonstrate SF, TF, and orient testing stim 

% movietype = 'sf';
% movietype = 'tf';
movietype = 'orient';



% defaults
pars.framerate = 60; % hz
pars.npeaks = 4; % number of white/black bars
pars.tempfreq = 1; % hz
pars.angle = 0; % degrees
pars.xlengthpix = 200; 
pars.ylengthpix = 200; 
pars.videolength = 5; % seconds
spacer_xwidth_pix = round(pars.xlengthpix / 5);
spacer_color = 0; % 0=black, 1=white
video_profile = 'Grayscale AVI';

pars.nframes_fullvideo = round(pars.videolength*pars.framerate);
spacerbar = repmat(spacer_color,pars.ylengthpix,spacer_xwidth_pix,pars.nframes_fullvideo);


switch movietype
    case 'sf'
        %% make sfreq variants
        outFileName = 'sf_grat_video';
        npeaks_list = [1 2 4 8];
        gratstack = [];
        for i = 1:length(npeaks_list)
            pars.npeaks = npeaks_list(i);
            if i ~= 1
                gratstack = [gratstack spacerbar];
            end
            new_stack = make_grat_stack(pars); 
            gratstack = [gratstack new_stack];
        end


    case 'tf'
        %% make tfreq variants
        outFileName = 'tf_grat_video';
        tfreq_list = [.5 1 2 4];
        gratstack = [];
        for i = 1:length(tfreq_list)
            pars.tempfreq = tfreq_list(i);
            if i ~= 1
                gratstack = [gratstack spacerbar];
            end
            new_stack = make_grat_stack(pars); 
            gratstack = [gratstack new_stack];
        end
    
    case 'orient'
        %% make angle variants
        outFileName = 'orient_grat_video';
        angle_list = [0 90 180 270];
        gratstack = [];
        for i = 1:length(angle_list)
            pars.angle = angle_list(i);
            if i ~= 1
                gratstack = [gratstack spacerbar];
            end
            new_stack = make_grat_stack(pars); 
            gratstack = [gratstack new_stack];
        end
        
    
end


%% write the video
myVideo = VideoWriter(outFileName,video_profile);
myVideo.FrameRate = framerate;
open(myVideo);
writeVideo(myVideo, gratstack);
close(myVideo);
clear gratstack pars new_stack myVideo






%%
function gratstack_out = make_grat_stack(pars)
    nframes_iter = round(pars.framerate / pars.tempfreq); % for 1 iteration
    repetitions = ceil(pars.nframes_fullvideo / nframes_iter); 
    
    if pars.angle == 0 % simplified for zero rotation
        gratstack_out = NaN(1,pars.xlengthpix,nframes_iter); % z = time
        for iframe = 1:nframes_iter
            gratstack_out(:,:,iframe) = sin([linspace(0,pars.npeaks*2*pi,pars.xlengthpix)] + iframe*2*pi/nframes_iter);
        end
        gratstack_out = repmat(gratstack_out,pars.ylengthpix,1,1); % expand into y
    else
        period_pix = round(pars.xlengthpix / pars.npeaks); 
        [y x]  = ndgrid(1:pars.ylengthpix, 1:pars.xlengthpix);
        th = deg2rad(pars.angle); 
        angle_grid = cos(th)*x + sin(th)*y;
        gratstack_out = NaN(pars.ylengthpix, pars.xlengthpix, nframes_iter);
        for iframe = 1:nframes_iter
            fraction_elapsed = iframe/nframes_iter; 
            gratstack_out(:,:,iframe) = sin( 2*pi * [angle_grid/period_pix+fraction_elapsed] ); 
        end
    end
    
    %%% fit to 0-1
    minstk = min(min(min(gratstack_out)));
    gratstack_out = gratstack_out - minstk;
    maxstk = max(max(max(gratstack_out)));
    gratstack_out = gratstack_out ./ maxstk;

    % loop the iteration
    gratstack_out = repmat(gratstack_out,1,1,repetitions);
    gratstack_out = gratstack_out(:,:,1:pars.nframes_fullvideo);


end