function RecCB(action, hfig)
if (nargin < 2)
        hfig = gcbf;
end
hobj = gcbo;
global gIsRecording;
ud = get(hfig,'UserData');
switch(action)
case 'UpdateInfo'
    axesType = get(findobj(hfig,'Tag','SourceMenu'),'Value');
    grid = [20 10 1235 600];
    [ax,width] = AxesLayout(hfig,grid,axesType);
    
    saveFile = get(findobj(hfig,'Tag','Save'),'Value');
    if(saveFile)
        set(findobj(hfig,'Tag','UsrHdr'),'Enable','on');
        set(findobj(hfig,'Tag','SaveAs'),'Enable','on');
        set(findobj(hfig,'Tag','File'),'Enable','on');
       
    else
        set(findobj(hfig,'Tag','UsrHdr'),'Enable','off');
        set(findobj(hfig,'Tag','SaveAs'),'Enable','off');
        set(findobj(hfig,'Tag','File'),'Enable','off');
    end
    
    stimulate = get(findobj(hfig,'Tag','Stimulate'),'Value');
    if(stimulate)
        set(findobj(hfig,'Tag','LoadStimFile'),'Enable','on');
        set(findobj(hfig,'Tag','RecordingMenu'),'Value',1);
        set(findobj(hfig,'Tag','StimulatorAxes'),'Visible','on');
    else
        set(findobj(hfig,'Tag','LoadStimFile'),'Enable','off');
        set(findobj(hfig,'Tag','RecordingMenu'),'Value',2);
        set(findobj(hfig,'Tag','StimulatorAxes'),'Visible','off');
    end
    
    recDuration = get(findobj(hfig,'Tag','RecordingMenu'),'Value');
    if (recDuration == 2)
        set(findobj(hfig,'Tag','RecordingTime'),'Enable','on');
    else
        set(findobj(hfig,'Tag','RecordingTime'),'Enable','off');
    end
    ud.ax = ax;
    ud.width = width;
    set(hfig,'UserData',ud);

    
case 'ChannelStatus'
    ud = get(hfig,'UserData');     
    ax = ud.ax;
    hobj = gcbo;
    set(findobj(hfig,'Tag','AllChannels'),'Value',0);
    indx = find([ax.axhandle] == hobj);
    if ((ax(indx).status)== 2)
       ax(indx).status = 1;
       set(ax(indx).labelhandle,'Color','k');
       ud.ax = ax;
       set(hfig,'UserData',ud);
       return
    end

    if ((ax(indx).status) == 1)
        ax(indx).status = 0;
        set(ax(indx).labelhandle,'Color','r');
        ud.ax = ax;
        set(hfig,'UserData',ud);
        return;
    end
    if ((ax(indx).status) == 0)
        ax(indx).status = 2;
        set(ax(indx).labelhandle,'Color','b');
        ud.ax = ax;
        set(hfig,'UserData',ud);
        return;
    end

case 'AllChannels'
    ud = get(hfig,'UserData');
    ax = ud.ax;
    allChannels = get(findobj(hfig,'Tag','AllChannels'),'Value');
    if (ax(1).status == 1 & allChannels)
        for i=1:prod(size(ax))
            ax(i).status = 1;
        end
    set([ax.labelhandle],'Color','k');
    ud.ax = ax;
    set(hfig,'UserData',ud);
    end 
    if (ax(1).status == 0 & allChannels)
        for i=1:prod(size(ax))
            ax(i).status = 0;
        end
    set([ax.labelhandle],'Color','r');
    ud.ax = ax;
    set(hfig,'UserData',ud);     
    end 
    if (ax(1).status == 2 & allChannels)
        for i=1:prod(size(ax))
            ax(i).status = 2;
        end
    set([ax.labelhandle],'Color','b');
    ud.ax = ax;
    set(hfig,'UserData',ud);
    end 
    
    
case 'SaveAsButton'
    [fname,pathname] = uiputfile('*.bin','Save data to file:');
                if (fname == 0)
                        disp('Operation cancelled');
                        return
                end
        fileName = [pathname,fname];
    hfile = findobj(hfig,'Tag','File');
    set(hfile,'String',sprintf('File:%s',fileName));
    ud.fileName = fileName;
    set(findobj(hfig,'Tag','Start'),'Enable','on'); % enable Stop button
    set(hfig,'UserData',ud);
    
    
case 'LoadStimFile'
    [fname,pathname]= uigetfile('*.mat','Load Stimulus File:');
    stimfileName = [pathname,fname];
    load(stimfileName);
    [xb,yb]= stairs(schedule(2,:),schedule(1,:));    
    hStimAxes = findobj(hfig,'Tag','StimulatorAxes');
    plot(xb,yb,'Parent',hStimAxes);
    %xlabel('Time (s)','Parent', hStimAxes);
    %ylabel('Valve Number','Parent',hStimAxes);
    lim=get(hStimAxes,'YLim');
    ud.hStimAxes = hStimAxes;
    ud.lim = lim;
    ud.schedule = schedule;
    ud.term = term;
    set(hfig,'UserData',ud);
    
    
case 'StartButton'
    saveFile = get(findobj(hfig,'Tag','Save'),'Value');
    stimulate = get(findobj(hfig,'Tag','Stimulate'),'Value');
    ud = get(hfig,'UserData'); 
    ax = ud.ax;
    pixWidth = ceil(ud.width);
    if(stimulate & ~isfield(ud,'hStimAxes'))
        errordlg('Must load a stimulus file','Input Error');
        return
    end
    if(saveFile & ~isfield(ud,'fileName'))
        errordlg('Must choose a file in order to save','Input Error');
        return
    end
    if(stimulate)
        hStimAxes = ud.hStimAxes;
        lim = ud.lim;
        schedule = ud.schedule;
        term = ud.term;
    else
        schedule = [];
    end
    [p.scanrate,recRange,recTime,dispRange,transferInterval,nblocks,p.usrheader] = RecCBGetData(hfig,schedule);
    if(recTime <= 0)
        errordlg('Must specify a recording duration','Input Error');
        return
    end
    
    gIsRecording = 1;
        set(findobj(hfig,'Tag','Stop'),'Enable','on'); % enable Stop button
    drawnow;
    dev = '/dev/comedi0';
    subdev = 0;
    activeIndex = find([ax.status] == 2);
    index = 1:length(ax);
    if (saveFile)
        index = find([ax.status] == 2 | [ax.status] == 1);
        fileName = ud.fileName;
        set(hobj,'Enable','off');        % disable Start button
    end
    disp('open')
    % Open the device
    if ~comedi_isopen(dev)
        comedi_open(dev);   
        end
        % Increase the buffer size to the maximum
        bufsz = comedi_get_max_buffer_size(dev,subdev);
        comedi_set_buffer_size(dev,subdev,bufsz);
    % Set up comedi command
    if (saveFile)
        cmd = make_comedi_cmd([ax(index).channel],p.scanrate,recRange);
    else
        cmd = make_comedi_cmd([ax.channel],p.scanrate,recRange);
    end    
        % Tweak command to comply with values actually supported by the board
        msg = 'fail';
        nfailed = 0;
        maxfail = 5;
        while (~strcmp(msg,'success') & nfailed < maxfail)
        [msg,cmd] = comedi_command_test(dev,cmd);  % Tunes parameters to realizable
        nfailed = nfailed + 1;
        end
        if (nfailed == maxfail)
        error('Command preparation did not succeed');
        end
        % Prepare the read buffer
    blockSize = p.scanrate * transferInterval;
    if (saveFile)
        readbuf = int16(zeros((length([ax(index).channel])),blockSize));
    else
        readbuf = int16(zeros((length([ax.channel])),blockSize));
    end
    voltageRange = comedi_get_range(dev,subdev,ax(index(1)).channel,recRange);
    p.voltageMax = voltageRange(2);
    p.voltageMin = voltageRange(1);
    displayRange = comedi_get_range(dev,subdev,ax(index(1)).channel,dispRange);
    % Calculate elements that are to be written to header
    p.scalemult = (p.voltageMax - p.voltageMin)/4095;
    p.scaleoff = p.voltageMin;
    
    if (saveFile)
        % Open the raw data file
        [fid, message] = fopen(fileName,'w');
        if (fid == -1)
            disp(message);
        end
        % Calculate elements that are to be written to header
        p.nscans = blockSize * nblocks;
        p.chan = [ax(index).channel];
        p.date = datestr(floor(now));
        p.time = datestr(rem(now,1));

        % Write the raw data header type
        positionAI = WriteAIHeader(fid,p,struct('AI',1));
    end
    
        nbufs = 0;
        curlength = 0;
    % Display time remaining of the read
    htimeLeft = findobj(hfig,'Tag','TimeLeftTxt');
    set(htimeLeft,'String',sprintf('Time Remaining: %2.0fs',recTime));
    drawnow
    timeLeft = recTime;
    if(stimulate)
        hMarker = line('XData',[0 0],'YData',lim,'Color','r','LineStyle','--','Parent',hStimAxes);
    end
    disp('read')
    if (comedi_command(dev,cmd) < 0)
        disp('Comedi: issuing command failed. Maybe a command is already running?');
        disp('You might fix this by calling comedi_cancel_and_flush for this device and then trying again.');
        error('');
        end
    if(stimulate)
        ScheduleStimulator(schedule);
    end
        while ((nbufs < nblocks) & gIsRecording)
        % Read raw data
        curlength = comedibuf_read(dev,readbuf,curlength);
        if (curlength < 0)
            error('Read command finished by itself!');
        end
        % Buffer is full
        if (curlength == prod(size(readbuf)))
            % Update and display time remaining of the read
            timeLeft = timeLeft - transferInterval;
            set(htimeLeft,'String',sprintf('Time Remaining: %2.1fs',timeLeft));
            drawnow;
            if(stimulate)
                delete(hMarker);
                hMarker = line('XData',[recTime-timeLeft recTime-timeLeft],'YData',lim,'Color','r','LineStyle','--','Parent',hStimAxes);
                drawnow;
            end
            % Save raw data to disk
            if (saveFile)
                fwrite(fid,readbuf(activeIndex,:),'int16');
            end 
            actualVoltage = p.scalemult * double(readbuf) + p.scaleoff;
            if(blockSize < 4*pixWidth)
                for i=1:prod(size([ax(index).channel]))
                    if (ax(index(i)).status == 2)
                        plot(actualVoltage(index(i),:),'Parent',ax(index(i)).axhandle,'Color','b');
                        set((ax(index(i)).axhandle),'YLim',displayRange);
                        set((ax(index(i)).axhandle),'YTick',[]);
                        set((ax(index(i)).axhandle),'XTick',[]);
                        drawnow;
                    end
                    if (ax(index(i)).status == 1)
                        plot(actualVoltage(index(i),:),'Parent',ax(index(i)).axhandle,'Color','k');
                        set((ax(index(i)).axhandle),'YLim',displayRange);
                        set((ax(index(i)).axhandle),'YTick',[]);
                        set((ax(index(i)).axhandle),'XTick',[]);
                        drawnow;
                    end
                end
            else
                envmm = envelopemem(actualVoltage,blockSize/pixWidth);
                xx = [1:pixWidth flipdim(1:pixWidth,2)];
                yy = [envmm.min flipdim(envmm.max,2)];
                for i=1:prod(size([ax(index).channel]))
                    if (ax(index(i)).status == 2)
                        fill(xx,yy,'b','EdgeColor','b','Parent',ax(index(i)).axhandle);
                        set((ax(index(i)).axhandle),'YLim',displayRange);
                        set((ax(index(i)).axhandle),'YTick',[]);
                        set((ax(index(i)).axhandle),'XTick',[]);
                        drawnow;
                    end
                    if (ax(index(i)).status == 1)
                        fill(xx,yy,'k','EdgeColor','k','Parent',ax(index(i)).axhandle);
                        set((ax(index(i)).axhandle),'YLim',displayRange);
                        set((ax(index(i)).axhandle),'YTick',[]);
                        set((ax(index(i)).axhandle),'XTick',[]);
                        drawnow;
                    end
                end
             end
            curlength = 0;
            nbufs = nbufs + 1;
        end  
    end
    
if(gIsRecording == 0 & stimulate)
    ScheduleStimulator(term);
end
 
% Close the data acquisition device
comedi_cancel_and_flush(dev,subdev);
comedi_close(dev);
disp('Done')


if (gIsRecording & stimulate)
  delete(hMarker);
end
    
% Close the file
    if (saveFile)
        status = fclose(fid);
        if (status < 0)
            error('File did not close');
        end   
    end   
    
case 'StopButton'
    disp('stop')
    gIsRecording = 0;
    set(hobj,'Enable','off');        % Disable Stop button
    drawnow;
end    

function [scanrate,recRange,recTime,dispRange,transferInterval,nblocks,usrheader] = RecCBGetData(hfig,schedule)
scanrate = str2num(get(findobj(hfig,'Tag','ScanRate'),'String'));
recRange = get(findobj(hfig,'Tag','RecordingRange'),'Value') - 1;
dispRange = get(findobj(hfig,'Tag','DisplayRange'),'Value') - 1;
transferInterval = str2num(get(findobj(hfig,'Tag','TransferInterval'),'String'));
saveFile = get(findobj(hfig,'Tag','Save'),'Value');
    if (saveFile)
        hUsrHdr = findobj(hfig,'Tag','UsrHdr');
        temphdr = get(hUsrHdr,'String')';
                ct = cellstr(temphdr');                % Pad header with returns at end of lines
                ret = double(sprintf('\n'));
                ctcat = strcat(ct,num2cell(char(ret*ones(length(ct),1))));
                usrheader = cat(2,ctcat{:})';        % Turn into a single column string
    else
        usrheader = [];
    end
recDuration = get(findobj(hfig,'Tag','RecordingMenu'),'Value');

    if (recDuration == 1)
        recTime = schedule(2,end);
        nblocks = recTime/transferInterval;
    elseif (recDuration == 2)
        recTime = str2num(get(findobj(hfig,'Tag','RecordingTime'),'String'));
        nblocks = recTime/transferInterval;
    else
        recTime = Inf;
        nblocks = Inf;
    end

return