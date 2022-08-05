function DoRWCB(action)
global gIsRecording
hfig = gcbf;
hobj = gcbo;
disp('enter')
switch (action)
case 'start'
    error
    disp('start')
        gIsRecording = 1;
        set(hobj,'Enable','off');        % disable Start button
        set(findobj(hfig,'Tag','StopRecButton'),'Enable','on'); % enable Stop button
%     drawnow
    p = get(gcbf,'UserData');
    dev = '/dev/comedi0';
    subdev = 0;
%     if ~comedi_isopen(dev)
%         comedi_open(dev);   % Open the device
%         end
%         % Increase the buffer size to the maximum
%         bufsz = comedi_get_max_buffer_size(dev,subdev);
%         comedi_set_buffer_size(dev,subdev,bufsz);
%         % Set up comedi command
%         cmd = make_comedi_cmd(p.chan,p.scanrate,p.gain);
        % Tweak command to comply with values actually supported by the board
        msg = 'fail';
        nfailed = 0;
        maxfail = 5;
   
%         while (~strcmp(msg,'success') & nfailed < maxfail)
%         [msg,cmd] = comedi_command_test(dev,cmd);  % Tunes parameters to realizable
%         nfailed = nfailed + 1;
%         end
%         if (nfailed == maxfail)
%         error('Command preparation did not succeed');
%         end
        % Prepare the read buffer
    block_size = 2 * p.nfreq * p.navg * p.nperblock;
        readbuf = int16(zeros(length(p.chan),block_size));
    % Opens the file
    if (p.savefile)
        [fid, message] = fopen(p.filename,'w');
        if (fid == -1)
            disp(message);
        end
    end
        % Now start reading
        nbufs = 0;
        curlength = 0;
    %Display time remaining of the read
    htime_left = findobj(hfig,'Tag','TimeLeftTxt');
    set(htime_left,'String',sprintf('Time Remaining: %2.0fs',p.tacq));
%     drawnow
    time_left = p.tacq;
        % Execute the command
        % After this, the board is actually running, and you need to empty the
        % read buffer sufficiently quickly
    disp('read')
%         if (comedi_command(dev,cmd) < 0)
%         disp('Comedi: issuing command failed. Maybe a command is already running?');
%         disp('You might fix this by calling comedi_cancel_and_flush for this device and then trying again.');
%         error('');
%         end
    zer = zeros(1,p.nfreq);
    xold = zer;
        while ((nbufs < p.nblocks) & gIsRecording)
%         curlength = comedibuf_read(dev,readbuf,curlength);
        curlength = curlength + ceil((prod(size(readbuf))-curlength)/2);
        if (curlength < 0)
            error('Read command finished by itself!');
        end
        % Buffer is full, do something with it
        if (curlength == prod(size(readbuf)))
            pause(1);
            %Update and display time remaining of the read
            htime_left = findobj(hfig,'Tag','TimeLeftTxt');
            time_left = time_left - p.tbuf;
            set(htime_left,'String',sprintf('Time Remaining: %2.1fs',time_left));
            % Save data to disk
            if (p.savefile)
                fwrite(fid, readbuf, 'int16');
            end
            % Compute the specgram
            [B,xold] = specgramwrap(readbuf,xold,2*p.nfreq);
            % Decimate sonogram in preparation for graphing
            %Bo = sngdecimate(abs(B),p.navg);
            %imagesc(Bo);
            
            
            
            
           
            
            
            curlength = 0;
            nbufs = nbufs + 1;
            drawnow
        end   
    end
    % Close the data acquisition
%     comedi_cancel_and_flush(dev,subdev);
%     comedi_close(dev);
    % Close the file
    if (p.savefile)
        [B1,xold] = specgramwrap(zer,xold,p.nfreq);
        status = fclose(fid);
        if (status < 0)
            error('File did not close');
        end   
    else
        set(hobj,'Enable','on');  % Enable start button
        set(findobj(hfig,'Tag','StopRecButton'),'Enable','off'); % Disable Stop button
    end
    disp('Done')
    
    %Plot data
%    if (p.savefile)
%          figure
%          fid = fopen(p.filename, 'r');
%          a = fread(fid, 'int16');
%          plot (a);
%          fclose(fid);
%          
%    end

case 'stop'
    disp('stop')
    gIsRecording = 0;
    set(hobj,'Enable','off');        % disable Stop button
        set(findobj(hfig,'Tag','StartRecButton'),'Enable','on'); % enable Start button
    drawnow
end