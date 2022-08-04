function htb = htbShow(filename, db, nEpoch)
%htbShow - Display an epoch from a database in an HTB file
%
% SYNOPSIS
%   htb = htbShow(filename, db, nEpoch);
%
% IN
%   filename            Name of HTB file
%   db                  Database number in file (=1, ...)
%   nEpoch              Epoch to display (1, ...)
%
% OUT
%   A 3d plot

fid = htbOpen(filename);            % Open HTB file
htb = htbGetHd(fid, db);            % Get HTB header
e = htbGetEp(htb, nEpoch);          % Get epoch data

esize = size(e);                    % =[rows columns]
if (esize(1) ~= 1)                  % Did we get any data?
    rotate3d on                     % Drag an axis for 3d rotation!
    t = sprintf('TEMPO HTB Database %d Epoch %d analyzed and plotted by Matlab.', db, nEpoch);
%   disp(t);                        % Display the title on Matlab's command line

%   mesh(e);                        % Plot as a wire diagram
%   colormap(zeros(64,3));          % Plot in black & white (for printing)

    surfc(e);                       % Plot as a surface with contour plot
%   surfl(e);                       % Plot as a surface
%   surf(e);                        % Plot as a surface
%   shading interp;                 % Plot as a smooth surface
    shading faceted;                % Plot with skeleton

    title(t);
    zlabel('Millivolts');
    xlabel('Channels');
    ylabel('Milliseconds');
    end

htbClose(fid);
return;