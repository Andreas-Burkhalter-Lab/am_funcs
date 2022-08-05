function DoRWCB(action)
global gIsRecording
switch(action)
case 'start'
        p = get(gcbf,'UserData');
        npix = 256;                        % This is currently not used in a meaningful way
                                                % (only used if graphing raw waveform in real time)

        % First, clear anything that might be in the axes
        %rwh = findobj(0,'Tag','SngWindow');
        %figure(rwh);
        htext = findobj(gcbf,'Tag','TimeLeftTxt');
        gIsRecording = 1;
        set(gcbo,'Enable','off');        % disable Start button
        set(findobj(gcbf,'Tag','StopRecButton'),'Enable','on'); % enable Stop button
        set(htext,'String','Recording...');
        imH = get(findobj(gcbf,'Tag','plotwindow'),'Children');
        %keyboard
        if (p.savefile)
                sng = RecordingWhis(p.chan,p.scanrate,p.nfreq,p.navg,p.nperblock,p.nblocks,p.gain,npix,imH,htext,p.video,p.usrheader,p.filename);
        else
                sng = RecordingWhis(p.chan,p.scanrate,p.nfreq,p.navg,p.nperblock,p.nblocks,p.gain,npix,imH,htext,p.video);
        end
        set(findobj(gcbf,'Tag','StopRecButton'),'Enable','off'); % disable Stop button
        gIsRecording = 0;
        sound(sin((1:1000)/4))                        % beep
        p.tacq = p.tacq*size(sng,2)/(p.nperblock*p.nblocks);        % correct the acq. time, if user hit stop
        if (p.savefile)                                        % Save the sonogram as a MAT file
                k = findstr(p.filename,'.bin');
                k1 = findstr(p.filename,filesep);
                if (isempty(k1))
                        k1 = 0;
                end
                if (length(k) ~= 1 | k < 2 | k1(end)+1 > length(p.filename))
                        name = '';
                else
                        name = [p.filename((k1(end)+1):k-1),'sng.mat'];
                end
                %[filename,pathname] = uiputfile(name);
                %pathname = pwd;
                pathname = p.filename(1:k1(end));
                filename = name;
                if (filename)
                        save([pathname,filename],'sng','p');
                end
        end
        figure
        ff = linspace(0,p.scanrate/2000,p.nfreq);
        fi = find(ff >= 40 & ff <= 90);
        mnsng = mean(sng(fi,:));
        mdmnsng = medfilt1(mnsng,11);
        plot(linspace(0,p.tacq,length(mnsng)),mdmnsng);
        axis tight;
        xlabel('Time (s)');
        ylabel('Mean power (40-90kHz)')
        if (p.savefile)
                title(p.filename);
        end
case 'stop'
        gIsRecording = 0;
        set(gcbo,'Enable','off');        % disable Stop button
case 'saveSng'
        % Note: this is currently not executed
        p = get(gcbf,'UserData');
        sng = getappdata(gcbf,'sng');
        %buttonname = questdlg('Save sonogram?','','Yes','No','Yes');
        buttonname = 'Yes';
        if strcmp('Yes',buttonname)
                k = findstr(p.filename,'.bin');
                k1 = findstr(p.filename,filesep);
                if (isempty(k1))
                        k1 = 0;
                end
                if (length(k) ~= 1 | k < 2 | k1(end)+1 > length(p.filename))
                        name = '';
                else
                        name = [p.filename((k1(end)+1):k-1),'sng.mat'];
                end
                [filename,pathname] = uiputfile(name);
                if (filename)
                        save([pathname,filename],'sng','p');
                end
        end
otherwise
        disp('DoRWCB: do nothing');
end
