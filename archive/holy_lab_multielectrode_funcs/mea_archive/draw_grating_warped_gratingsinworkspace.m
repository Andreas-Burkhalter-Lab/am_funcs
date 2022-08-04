% Present grating warped to resemble a sinusoidal grating on a sphere
% projected onto the 2D screen, as per Marshel et al. 2011
%%% This version creates the gratings as variables in the matlab
%%% workspace, which overloads memory if there are too many trials (more
%%% than ~20). For larger numbers of trials, use newer version that loads
%%% from saved gratings files. 
%%% last updated 6/19/15 on msi

%%% prerender reduces ifi to 1/60hz, but requires long time to prerender

clear
%% Setup
%%% vector elements specify the outer-grating parameter for a different iteration;   
%%% inner grating sf and tf are fixed; inner and outer amp is fixed;
%%% gratings always iso-oriented
prerender = 1; % make textures before entering the Flip loop; if off, only take parameters from the first element
z0 = 300;  % length in pix of perpendicular from eye to screen (? =zero point of x and y)... need to enable declaration as inches
inner.sf = 200;
inner.tf = 0.5;
inner.amp = 1;
inner.diam = 30;
outer.sf = [100 500 200 800]% 50 300 600]; % cycles per degree; peak ~0.05 in Gao..... currently not cycles per degree
outer.tf = [10 5]% 0.2 1 3 0.8 0.5]; % hz
outer.amp = 1;
outer.diam = [35 50 70 80]% 90 100 110];
theta = [0 60]% 120 180 240 300];   %%%% degrees; positive is clockwise
dur = 2; % seconds
ifi = 1/60;
isi = 0.5; % interstimulus interval in seconds
nframes = round(dur/ifi);
stim_center = [800 600]; % projective center for drawing aperture; not the 2D stim center in pixels
% stim_center = follow_cursor;
diam = [13 20 30 50 70 100]; %stim diameter in degrees

Screen('CloseAll'); %% clear all previously drawn textures to free up memory
myScreen = max(Screen('Screens'));
[win , winrect] = Screen(myScreen,'OpenWindow',BlackIndex(myScreen));
screenstats = Screen('Resolution',win);
scrnperpx = round(screenstats.width/2); % perpendicular to screen through eye must pass through this x-value; usually screen center
scrnperpy = round(screenstats.height/2); % perpendicular to screen through eye must pass through this y-value; usually screen center
[xcentermesh ycentermesh] = meshgrid(1:screenstats.width,1:screenstats.height);
perp2stimcenterx = scrnperpx - stim_center(1);
perp2stimcentery = scrnperpy - stim_center(2);
perp2meshx = scrnperpx-xcentermesh;
perp2meshy = scrnperpy-ycentermesh;
[xstraight ystraight] = meshgrid(-screenstats.width/2 : screenstats.width/2, -screenstats.height/2 : screenstats.height/2);

% Construct stimuli for each frame.
tic
for patch = [1 2]   % for the center grating and annulus
            %%%% if need to make more efficient, just compute the grating in
        %%%% the relevant rect..... but this won't help when the grating is
        %%%% large anyway
    if prerender
        switch patch 
            case 1
                % Make a table for listing inner-grating texture handles. 
                rep_tex = cell(length(theta),1); 
                tex_inner = table(rep_tex, theta','VariableNames',{'tex','Angle'});
                sfrad = deg2rad(inner.sf);
                
                % Make the inner-grating aperture.
                clear aperture
                apt_fullscreen = deg2rad(inner.diam/2) > acos( (z0^2 + perp2stimcenterx*perp2meshx + perp2stimcentery*perp2meshy) ./...
                    sqrt( (z0^2+perp2stimcenterx^2+perp2stimcentery^2)*(z0^2+perp2meshx.^2+perp2meshy.^2)) );
                [apty aptx] = find(apt_fullscreen); % coordinates of pixels covered by the aperture
                inner.rect_aperture = [min(aptx) min(apty) max(aptx) max(apty)];
                aperture(:,:,4) = 255*double(apt_fullscreen(min(apty):max(apty), min(aptx):max(aptx))); % draw only the nonzero rect of aperture into alpha; RGB=0
                inner.apertex = Screen('MakeTexture', win, aperture);
                
                % Make the inner gratings. Gratings move along the yrotate dimension.
                for Angle = 1:length(theta)  % could calculate these values just once for each angle, then call them for inner and outer
                    thetarad = deg2rad(theta(Angle));
                    xrotate = sin(-thetarad)*ystraight + cos(-thetarad)*xstraight;
                    yrotate = cos(-thetarad)*ystraight + sin(-thetarad)*xstraight;
                    eye2screen = pi/2 - acos(yrotate ./ sqrt(z0^2 + xrotate.^2 + yrotate.^2)); % angles in radians; pi/2 - acos is faster than asin
                    for frame = 1:nframes
                        grating = 0.5*WhiteIndex(win) + 0.5*WhiteIndex(win)*inner.amp*cos(2*pi*sfrad*(eye2screen) - frame*ifi*inner.tf); % from Marshel et al.
                        match = tex_inner.Angle == theta(Angle);
                        tex_inner.tex{match}(frame) = Screen('MakeTexture', win, grating);
                    end
                end
            case 2
                % Make a table for outer-grating listing texture handles. Create parameter values iteratively, rather than as a single vector above.
                totaltrials = length(theta)*length(outer.diam)*length(outer.sf)*length(outer.tf);
                nans = NaN(totaltrials,1);
% % % %                 tex_outer = table(cell(size(nans)),nans,nans,nans,nans,'VariableNames',{'tex','diam','Angle','sf','tf'});
                tex_outer = table(cell(size(nans)),nans,nans,nans,nans,'VariableNames',{'grating','diam','Angle','sf','tf'});
                outer.apertex = NaN(1,length(diam)); % initialize
                
                count = 0; % to index into successive tex_outer rows
                for diam = 1:length(outer.diam)
                    % Make the outer-grating aperture.
                    clear aperture
                    apt_fullscreen = deg2rad(outer.diam(diam)/2) > acos( (z0^2 + perp2stimcenterx*perp2meshx + perp2stimcentery*perp2meshy) ./...
                        sqrt( (z0^2+perp2stimcenterx^2+perp2stimcentery^2)*(z0^2+perp2meshx.^2+perp2meshy.^2)) );
                    [apty aptx] = find(apt_fullscreen); % coordinates of pixels covered by the aperture
                    outer.rect_aperture{diam} = [min(aptx) min(apty) max(aptx) max(apty)];
                    aperture(:,:,4) = 255*double(apt_fullscreen(min(apty):max(apty), min(aptx):max(aptx))); % draw only the nonzero rect of aperture into alpha; RGB=0
                    outer.apertex(diam) = Screen('MakeTexture', win, aperture);
                    
                    % Make the outer gratings. Gratings move along the yrotate dimension.
                    for Angle = 1:length(theta) % could calculate these values just once for each angle, then call them for inner and outer
                        thetarad = deg2rad(theta(Angle));
                        xrotate = sin(-thetarad)*ystraight + cos(-thetarad)*xstraight;
                        yrotate = cos(-thetarad)*ystraight + sin(-thetarad)*xstraight;
                        eye2screen = pi/2 - acos(yrotate ./ sqrt(z0^2 + xrotate.^2 + yrotate.^2)); % angles in radians; pi/2 - acos is faster than asin
                        for sf = 1:length(outer.sf)
                            sfrad = deg2rad(outer.sf(sf)); 
                            for tf = 1:length(outer.tf)
                            count = count+1;
                                for frame = 1:nframes
                                    tex_outer.diam(count) = outer.diam(diam);
                                    tex_outer.Angle(count) = theta(Angle);
                                    tex_outer.sf(count) = outer.sf(sf);
                                    tex_outer.tf(count) = outer.tf(tf);
                                    grating = 0.5*WhiteIndex(win) + 0.5*WhiteIndex(win)*outer.amp*...
                                        cos(2*pi*sfrad*(eye2screen) - frame*ifi*outer.tf(tf)); % from Marshel et al.
% % % %                                     tex_outer.tex{count}(frame) = Screen('MakeTexture', win, grating);
                                    tex_outer.grating{count}{frame} = grating; % takes too much memory
                                end
                            end
                        end
                    end
                end
                tex_outer = tex_outer(randperm(height(tex_outer)),:); %% shuffle the trial order
        end
    end
end
                                  
toc
%% Present stimuli
for trial = 1:height(tex_outer)
    for frame = 1:nframes
        for patch = [2 1] % Draw the annulus first so that the inner grating covers it. 
            if prerender   %% maybe add functionality for ~prerender
                if patch==2;            
% % % %                     gratingtex = tex_outer.tex{trial}(frame);
                    grating = tex_outer.grating{trial}(frame);
                    gratingtex = Screen('MakeTexture', win, grating);
                    match = outer.diam==tex_outer.diam(trial);
                    apertex = outer.apertex(match); % get the aperture of the appropriate diameter
                    rect_aperture = outer.rect_aperture{match};
                    Angle = tex_outer.Angle(trial);
                else
                    match = tex_inner.Angle == Angle;
                    gratingtex = tex_inner.tex{match}(frame);  %% match the inner grating angle to the outer grating angle
                    apertex = inner.apertex;
                    rect_aperture = inner.rect_aperture;
                end
            end

            % Disable alpha-blending, restrict following drawing to alpha channel:    
            Screen('Blendfunction', win, GL_ONE, GL_ZERO, [0 0 0 1]);

            % Clear 'dstRect' region of framebuffers alpha channel to zero:
            Screen('FillRect', win, [0 0 0 0], rect_aperture); % maybe take smaller of winrect and this rect
            % Fill the region specified by apertex with an alpha value of 255:
            Screen('DrawTexture', win, apertex, [], rect_aperture); 

            % Enable DeSTination alpha blending and reenable drawing to all
            % color channels. Following drawing commands will only draw there
            % the alpha value in the framebuffer is greater than zero, ie., in
            % our case, inside the circular 'dstRect' aperture where alpha has
            % been set to 255 by our 'FillOval' command:
            Screen('Blendfunction', win, GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA, [1 1 1 1]);

            Screen('DrawTexture', win, gratingtex, rect_aperture, rect_aperture) %%    % draw the warped grating
        end
        vbl(frame) = Screen('Flip',win);
    end
    Screen('FillRect',win,BlackIndex(win)); Screen('Flip',win);
    pause(isi);
end




