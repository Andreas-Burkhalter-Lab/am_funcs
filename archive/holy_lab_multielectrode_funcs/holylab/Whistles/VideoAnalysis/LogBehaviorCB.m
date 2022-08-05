function LogBehaviorCB(action,hfig)
if (nargin < 2)
        hfig = gcbf;
end
switch(action)
case 'Key'
        sound(sin((1:50)/4))                        % beep
        c = get(hfig,'CurrentCharacter');
        [status,tc] = SonySerialGate(4);        % Read current time code
        if ~status
                return;                                        % if the read is bad, quit here
        end
        t = getappdata(hfig,'t');
        tcur = tc(4)*3600+tc(3)*60+tc(2)+tc(1)/30 - t.start;
        k = getappdata(hfig,'BehavKeys');
        nbehavs = length(k);
        if (strcmp(c,' '))
                % It's a plot-forcing key (space bar)
                toggle = getappdata(hfig,'toggle');
                tindx = find(toggle(2,:) < tcur);        % Keep only past values (permit editing)
                toggle = toggle(:,tindx);
                setappdata(hfig,'toggle',toggle);
                val = ToggleToStair(toggle,tcur,nbehavs);
                PlotResp(val);
                return;
        end
        match = strcmp(c,k);
        indx = find(match);
        if (length(indx) ~= 1)
                % Check to see if it's a camcorder control button
                k = getappdata(hfig,'VideoKeys');
                indx = findstr(c,k);
                if (length(indx) ~= 1)
                        return
                end;
                vidcmd = getappdata(hfig,'VideoCommand');
                status = SonySerialGate(vidcmd(indx));
                return;
        end
        % It is a behavior-logging key. Update records.
        toggle = getappdata(hfig,'toggle');
        tindx = find(toggle(2,:) < tcur);        % Keep only past values (permit editing)
        toggle = toggle(:,tindx);
        toggle(:,end+1) = [indx;tcur];
        setappdata(hfig,'toggle',toggle);
        val = ToggleToStair(toggle,tcur,nbehavs);
        PlotResp(val);
case 'Plot'
        toggle = getappdata(hfig,'toggle');
        t = getappdata(hfig,'t');
        nbehavs = length(getappdata(hfig,'BehavKeys'));
        val = ToggleToStair(toggle,t.cur,nbehavs);
        PlotResp(val);
end

function PlotResp(val)
n = length(val);
hline = findobj(gca,'Tag','BehavLine');
hdots = findobj(gca,'Tag','BehavDots');
offset = -0.25;
amplitude = 0.5;
ye = zeros(1,n);
xe = ye;
for i = 1:n
        [xx,yy] = stairs(val{i}(2,:),val{i}(1,:)*amplitude+offset + i);
        if (isempty(hline))
                line(xx,yy,'Tag','BehavLine','EraseMode','background');
        else
                set(hline(i),'XData',xx,'YData',yy);
        end
        xe(i) = xx(end);
        ye(i) = yy(end);
end
if (isempty(hdots))
        hdots = line(xe,ye,'Tag','BehavDots','EraseMode','background',...
                'LineStyle','none','Marker','.','MarkerSize',9,'Color','r');
else
        set(hdots,'XData',xe,'YData',ye);
end        
drawnow

