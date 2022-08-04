%% Timer callback function that measures mouse displacement between time intervals
% Uses and modifies variables within the scope of draw_grating_warped. 
%%% last edited 9/30/15
function trackmouse_callback(h)
    global mouse_movements mousevars   % takes about 0.3ms
    mousevars.mouseScanIndex = mousevars.mouseScanIndex + 1;
    mousevars.scantime(mousevars.mouseScanIndex) = now; % time of this scan for syncing later
    [xNow yNow] = GetMouse;
    
    
    %%% below commands were used when mouse tracking was controlled in the same
    %%% Matlab window as stimulus presentation; now using 2 Matlab windows,
    %%% syncing via 'now'
    %     mousevars.mouseTic = tic;
    %     mousevars.timeSincePreviousMouseScan(mousevars.mouseScanIndex) = toc(mousevars.mouseTic);
    
    % If the mouse moved since the last scan, record this scan number (col 1)
    % within the trial (col 1) and the x (col 2) and y (col 3) displacements. 
    if any([xNow yNow] ~= mousevars.mousePre); % if there was a mouse displacement
        mouse_movements(mousevars.mouseScanIndex,1:2) = [xNow yNow]-mousevars.mousePre;
    end
    
    % Recenter mouse position if necessary and set current mouse position
    % as mousePre for the next scan. 
    if abs(xNow - mousevars.mouseScreenCenter(1)) > mousevars.xResetDist ||...
            abs(yNow - mousevars.mouseScreenCenter(2)) > mousevars.yResetDist
        SetMouse(mousevars.mouseScreenCenter(1),mousevars.mouseScreenCenter(2));
        mousevars.mousePre = mousevars.mouseScreenCenter;
    else
        mousevars.mousePre = [xNow yNow];
    end
    
end