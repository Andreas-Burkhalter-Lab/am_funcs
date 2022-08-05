function imstimplay(ip,options)
% IMSTIMPLAY: play stimulus-triggered movie snippets
% Syntax:
%  imstimplay(ip,options)
% where
%  ip is an imphys structure;
%  options is a structure with the following fields:
%    flush (default 0): the number of the flush valve;
%    trange (default [-12 20]): the time, in seconds, to include around a
%      valve transition;
%    valvenum (default all): valves to analyze;
%    datatypebg (defaults to datatype of ip): a string specifying the
%      datatype of the background frames (e.g., 'int16');
%    datatyperesp (defaults to datatype of ip): a string specifying the
%      datatype of the response frames;

% Copyright 2005 by Timothy E. Holy & Jason Guo

% Notes: the datatype stuff is all messed up. With spatial filtering, one
% can't use integer types, since sub-integer values are very meaningful.
% Single seems to be the best choice.

% Parse options
  if (nargin < 2)
    options = struct;
  end
  options = ispoptions(ip,options);
  if isfield(options,'valvenum')
    valvenum = options.valvenum;
  else
    valvenum = setdiff([ip.stimulus],options.flush);
  end
  % Find valve changes and collect ranges for background & responses
  [ipindx,identity] = imintervalsfromstim(ip,options.trange,valvenum);
    
  % Show an image to get crop region
  hCropRect = imselrect(imphysfetch(ip(1)));
  
  % Set up listdlg and do the main loop
  cellidentity = num2cell(identity);
  for i = 1:length(cellidentity)
    stridentity{i} = num2str(cellidentity{i});
  end
  ok = 1;
  choiceindx=1;
  while ok
    [choiceindx,ok] = listdlg('PromptString','Select a valve transition:',...
                              'SelectionMode','multiple',...
                              'InitialValue', choiceindx, ...
                              'ListString',stridentity);
    if ok
      if(isfield(options, 'implayer') && ishandle(options.implayer) )
         % close(options.implayer);
         filter_hsize=getappdata(options.implayer, 'filter_hsize');
         filter_sigma=getappdata(options.implayer, 'filter_sigma');
      else
         filter_hsize=25;
         filter_sigma=8;
      end
      % Create the spatial filter
      smoothfilt = fspecial('gaussian',filter_hsize, filter_sigma); 
      
      ipsub = cell(1,length(choiceindx));
      ipraw = ipsub;
      nSelectedValves=length(choiceindx);
      for i = 1:length(choiceindx)
        % Pick the period that we'll be working with
        ipcur = ip(ipindx{choiceindx(i)});
        % Retrieve crop rectangle (using a function for this might allow
        % more flexible notions of cropping, e.g., when we are using a
        % stack rather than single frames)
        ipcur = imcroprect(ipcur,hCropRect);
        % The background is computed from the initial period with
        % flush. Determine the period with flush, making sure that we
        % don't include stacks which occur after the stimulus valve turns
        % off again.
        bgindx = find([ipcur.stimulus] == options.flush);
        jumpindx = find(diff(bgindx) > 1);   % Interrupted by stimulus?
        if ~isempty(jumpindx)
          % bgindx = bgindx(1:jumpindx(1)+1);  % Just the initial period
          bgindx = bgindx(1:jumpindx(1) );  % Just the initial period
        end
        if isempty(bgindx)
          error('No background stacks available');
          % Or do we want to continue with no background subtraction?
        end
        % Compute the background
        tmp = double(imphysfetch(ipcur(bgindx(1))));
        for j = 2:length(bgindx)
          tmp = tmp + double( imphysfetch(ipcur(bgindx(j))) );
        end
        background = single(tmp / length(bgindx));   %TEH this used to have datatypebg
        % Now compute the background-subtracted stimulus stacks
        % stimindx = setdiff(1:length(ipcur),bgindx);
        stimindx = 1:length(ipcur);
        ipsub{i} = imphyscopy(ipcur(stimindx),{'dimension','timing','info', 'storage'});
        ipraw{i} = ipsub{i};
        % Might we need a progress bar here?
        nimages=length(ipsub{i});
        for j = 1:nimages
          imtmp = imphysfetch(ipcur(stimindx(j)));
          ipraw{i}(j).image = feval(options.datatyperesp,imtmp) ;
          ipraw{i}(j).background = background; % sacrafy some mem for flexibility of implayer
          ipsub{i}(j).image = ipraw{i}(j).image;
          % ipsub{i}(j).image = ...
          %   imfilter(single(imtmp) - background,smoothfilt);
          %ipsub{i}(j).image = feval(options.datatyperesp,single(imtmp) - background);
          
          [tt,tt_basename,tt_extname] = fileparts(ipcur(stimindx(j)).imfile);
          if(mod(j,5)==1 || j==nimages)
             progress_bar(struct('progress', j+(i-1)*nimages, ...
                                 'max', nimages*nSelectedValves, ...
                                 'what', ['loading ' [tt_basename,tt_extname] ' ...']));
          end

        end
        % todo: verify if the ipsub(i) has same stimulus?
      end
      
      % Create the player window 
      % This would basically be imanalyze, rewritten to accept imphys
      % structure as the input.
      % I might also recommend a couple of UI changes, esp. with regards
      % to the contrast settings, but we can get to that later.
      
      if(isfield(options, 'implayer') && ishandle(options.implayer) )
         % close(options.implayer);
         set(0, 'currentFigure', options.implayer);
         implayer('update_ip', options.implayer, [], [], ...
           cat(2,ipsub{:}), cat(2,ipraw{:}));
      else
         hImplayer=implayer(cat(2,ipsub{:}), cat(2,ipraw{:}));
         options.implayer=hImplayer;
         setappdata(options.implayer, 'filter_hsize', filter_hsize);
         setappdata(options.implayer, 'filter_sigma', filter_sigma);
      end
      set(options.implayer, 'waitstatus', 1); % 1: waitstatus is waiting
      waitfor(options.implayer, 'waitstatus', 0); % 0: not waiting
    end
  end
  
  
function options = ispoptions(ip,options)
  if ~isfield(options,'flush')
    options.flush = 0;
  end
  if ~isfield(options,'trange')
    options.trange = [-12 20];
  end
  % todo:
  options.datatypebg='int16';
  options.datatyperesp='int16';
%   if ~isfield(options,'datatypebg')
%     options.datatypebg = ip(1).datatype;
%   end
%   if ~isfield(options,'datatyperesp')
%     options.datatyperesp = ip(1).datatype;
%   end
%   
