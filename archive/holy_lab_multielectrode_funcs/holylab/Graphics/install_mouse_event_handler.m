function install_mouse_event_handler(obj, event_type, funcHandle)
% INSTALL_MOUSE_EVENT_HANDLER: install a mouse event handler to an object
% syntax:
%    install_mouse_event_handler(obj, event_type, funcHandle)
% pre: 
%    obj: the window to send event
%    event_type: one of 'down', 'up', 'move'/'motion'
%               (actually, only the first letter matters, and is case insensitive)
%    funcHandle: a string of matlab expression or a function handle 
%                that denotes which func to call when event occurs
%         its prototype is: 
%            handled=f(sender, event_args)
%         where
%             sender:     the object that receives the event;
%             event_args: a struct who has one field (event_type);
%             handled:    1 if you don't want the event to be
%                         handled by a handler installed later.
%         if you don't want parent window handles the event, install
%            handler for children ealier than for parent and
%            return 0 in children's event handler.
% eg:
%    figure; plot(30:40);
%    install_mouse_event_handler(gca, 'down', 'disp(''mouse down''); disp(get(gca, ''currentpoint''));') 
%    install_mouse_event_handler(gca, 'up', 'disp(''mouse up''); disp(get(gca, ''currentpoint''));')
%    uicontrol('Style', 'pushbutton', 'String', 'test','Position', [20 150 100 70]);
%    install_mouse_event_handler(gca, 'move', @test_event_data)
% 
% Note: if you need to install handlers on many objects, this
% mouse_event_handler code comes with some performance penalties. In such
% cases it may be better to use the basic Matlab functionality.

   hFig=get_parent_fig(obj);
   
   % to make is_mouse_over() work, check the units of position
   if(~has_good_unit(hFig) || ~has_good_unit(obj))
      % todo: should return a true/false instead of error
      error('unsupported unit setting in parent figure or the object');
   end

   % normalize event_type:
   switch(upper(event_type(1))) % only the first letter matters
      case 'D'
         event_type='down';
      case 'U'
         event_type='up';
      case 'M'
         event_type='motion';
      otherwise
         error(['unrecognized event_type: ' event_type]);
   end % switch

   % install the handlers to parent fig if necessary,
   % down and up must be installed; and motion is installed only when event_type is motion
   all_event_types={'down', 'up'};
   if(strcmp(event_type, 'motion')) all_event_types(end+1)={event_type}; end
   
   for idx=1:length(all_event_types)
      t_event_type=all_event_types{idx};
      handler_old=get(hFig, ['WindowButton' upper(t_event_type(1)) t_event_type(2:end) 'Fcn']);
      % since the eventdata passed by matlab is empty, we have to
      % use 3 handlers to distinguish 3 event types
      handler_new_str=['@mouse_event_dispatcher_' t_event_type];
      handler_new=eval(handler_new_str); % to use subfunctions in a file as handler, 
                                         % we must use a function handler.
                                         % OR: use str2func
      if(isempty(handler_old))
         set(hFig, ['WindowButton' upper(t_event_type(1)) t_event_type(2:end) 'Fcn'], handler_new);
      else
         if(~(isa(handler_old, 'function_handle') && isequal(handler_old,handler_new) ) )
            error(['the parent figure of obj ' num2str(obj) ...
                   ' has a handler for mouse ' t_event_type]);
         end
      end
   end
   
   % save the info of specific handler to obj in the parent figure's appdata 
   % w/ name mouse_event_handler_g :
   % note: use mouse_event_handler_g to save such info:
   %     mouse_event_handler_g.down is a 2xn cell array for mouse_down;
   %     mouse_event_handler_g.down(1,:) is objects;
   %     mouse_event_handler_g.down(2,:) is handlers;
   
   mouse_event_handler_g=getappdata(hFig, 'mouse_event_handler_g'); 
   % use cell array since func handle cann't be put into a matrix
   all_event_types={'down', 'up', 'motion'};
   for idx=1:length(all_event_types)
      t_event_type=all_event_types{idx};
      if(~isfield(mouse_event_handler_g, t_event_type))
	 mouse_event_handler_g.(t_event_type)={}; % to make "(:,end+1)=" works, and make 
						  % mouse_event_dispatcher() case 'd' 
						  % efficient
      end
   end
   
   if(~isfield(mouse_event_handler_g, 'is_mouse_down'))
      mouse_event_handler_g.is_mouse_down=0;
   end
   if(~isfield(mouse_event_handler_g, 'obj_when_mouse_down'))
      mouse_event_handler_g.obj_when_mouse_down=[];
   end
   if(~isfield(mouse_event_handler_g, 'pos_when_mouse_down'))
      mouse_event_handler_g.pos_when_mouse_down=[];
   end
   
   mouse_event_handler_g.(event_type)(:,end+1)={obj, funcHandle};

   setappdata(hFig, 'mouse_event_handler_g', mouse_event_handler_g); % save back the data shared 
								     % within the figure
   
   
% since the eventdata passed by matlab is empty, we have to
% use 3 handlers to distinguish 3 event types   
function mouse_event_dispatcher_down(sender, eventdata)
% the first two args are required by matlab   
   mouse_event_dispatcher(sender, eventdata, 'down');
   
function mouse_event_dispatcher_up(sender, eventdata)
   mouse_event_dispatcher(sender, eventdata, 'up');

function mouse_event_dispatcher_motion(sender, eventdata)
   mouse_event_dispatcher(sender, eventdata, 'motion');
   
% this is the one who really dispatches events to obj
function mouse_event_dispatcher(sender, eventdata, event_type)
% note: even user didn't install any handler for mouse down, we need update
%       is_mouse_down/pos_when_mouse_down/obj_when_mouse_down so that 
%       it still works when user only installed mouse_up handler.
% 
   hFig=sender; % the sender is the figure, when we dispatch 
		% the event, we need pass the "real sender"
   mouse_event_handler_g=getappdata(hFig, 'mouse_event_handler_g'); 
   
   switch(event_type(1)) % only the first letter matters
      case 'd'
         mouse_event_handler_g.is_mouse_down=1;
         mouse_event_handler_g.pos_when_mouse_down=get(0, 'PointerLocation');
         mouse_event_handler_g.obj_when_mouse_down=[]; % will be filled later
         
         % put all objects that installed mouse events into the variable "objects":
         objects=[];
         if(~isempty(mouse_event_handler_g.down))
            objects=[objects cell2mat(mouse_event_handler_g.down(1,:))];
         end
         if(~isempty(mouse_event_handler_g.up))
            objects=[objects cell2mat(mouse_event_handler_g.up(1,:))];
         end
         if(~isempty(mouse_event_handler_g.motion))
            objects=[objects cell2mat(mouse_event_handler_g.motion(1,:))];
         end
         objects=unique(objects);
         
         for idx=1:length(objects)
            if(is_mouse_over(objects(idx)))
               mouse_event_handler_g.obj_when_mouse_down(end+1)=objects(idx);
            end
         end
         setappdata(hFig, 'mouse_event_handler_g', mouse_event_handler_g); % save back
      case 'u'
         mouse_event_handler_g.is_mouse_down=0;
         setappdata(hFig, 'mouse_event_handler_g', mouse_event_handler_g); % save back
   end % switch
   tInfo=mouse_event_handler_g.(event_type);
   
   for idx=1:size(tInfo,2)
      tObj=tInfo{1,idx};
      if(is_event_receiver(tObj, event_type, mouse_event_handler_g)==1) % 
         tFuncHandle=tInfo{2,idx};
         if(isstr(tFuncHandle))
            % if the hanlder is a string of matlab expression
            eval(tFuncHandle);
         else
            tParam=struct('event_type', event_type);
            if(tFuncHandle(tObj, tParam)) % 1st arg: sender; 2nd arg: event args
               break;
            end % if, handled, don't look for next handler
         end % else, the handler is a function handle
      end % if, mouse is over the obj
   end % for, each handler

   
function result=is_event_receiver(obj, event_type, aMouse_event_handler_g)
   switch(event_type(1)) % only the first letter matters and is case sensitive
      case 'd'
         result=is_mouse_over(obj);
      case 'u'
         oldPos=aMouse_event_handler_g.pos_when_mouse_down;
         newPos=get(0, 'PointerLocation');
         if isequal(oldPos,newPos) % if not moved
            result=is_mouse_over(obj);
         else
            tt=find(aMouse_event_handler_g.obj_when_mouse_down==obj);
            if(isempty(tt))
               result=0;
            else
               result=1;
            end
         end % else, mouse position moved
      case 'm'
         if(aMouse_event_handler_g.is_mouse_down)
            tt=find(aMouse_event_handler_g.obj_when_mouse_down==obj);
            if(isempty(tt))
               result=0;
            else
               result=1;
            end
         else
            result=is_mouse_over(obj);
         end % else, mouse button is not down
      otherwise
         error('error in code');
   end % switch

