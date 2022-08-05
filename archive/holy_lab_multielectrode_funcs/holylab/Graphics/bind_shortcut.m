function bind_shortcut(obj, shortcut, funcHandle)
% BIND_SHORTCUT(obj, shortcut, funcHandle)
%  obj: a figure or uicontrol handle; shortcut: case insensitive, use + to separate modifier and key
%  funcHandle: 
%     a str, cell array or a function handle of prototype
%         f(sender, event_args)
%      OR f(sender, event_args, arg1, arg2, ...)
%     where: sender: the object; event_args: current is [].
% EG: figure; bind_shortcut(gcf, 'ctrl+b', 'disp(''key ctrl+b pressed'')')
%     figure; bind_shortcut(gcf, 'f1', {@show_help, 'some title', 'some msg'});
%
% SEE: keyPressFcn property.
% SEE: test_keypress()

   type=get(obj, 'Type');
   if(~(strcmp(type, 'figure') || strcmp(type, 'uicontrol')))
      error(['object ' num2str(obj) 'cannot be binded to shortcut']);
   end
   
   handler_old=get(obj, 'keyPressFcn');
   handler_new=str2func('shortcut_dispatcher');
   if(isempty(handler_old))
      set(obj, 'keyPressFcn', handler_new);
   else
      if(~(isa(handler_old, 'function_handle') && isequal(handler_old,handler_new) ) )
	 error(['the object ' num2str(obj) ...
		' has a handler for key press']);
      end
   end
   
   shortcut=lower(shortcut);
   keys=cellstr(split_str(shortcut, '+'));
   for idx=1:length(keys)
      switch(keys{idx})
         case {'c','ctrl'}
            keys{idx}='control';
      end
   end
   
   if(shortcut(end)=='+')
      bind.key='+';
   else
      bind.key=keys{end};
   end
   if(length(keys)<=1)
      bind.mod=cell(1,0);
   else
      bind.mod=sort(keys(1:end-1));
   end
   bind.action=funcHandle;
   
   shortcut_dispatcher_g=getappdata(obj, 'shortcut_dispatcher_g');
   shortcut_dispatcher_g{end+1}=bind;
   setappdata(obj, 'shortcut_dispatcher_g', shortcut_dispatcher_g);
   
   
function shortcut_dispatcher(sender, eventdata)
% TODO: may speed it up by hashing shortcut
   obj=sender;
   shortcut_dispatcher_g=getappdata(obj, 'shortcut_dispatcher_g');
   mod=sort(eventdata.Modifier);
   key=eventdata.Key;
   % if(isequal(key, 'b')) keyboard;  end
   for idx=1:length(shortcut_dispatcher_g)
      bind=shortcut_dispatcher_g{idx};
      if(isequal(bind.key, key) && isequal(bind.mod, mod))
	 action=bind.action;
	 if(isstr(action))
	    eval(action);
         else
            if(iscell(action))
               tfunc=action{1};
               tfunc(obj, [], action{2:end});
            else
               action(obj, []);
            end
	 end
      end
   end

   
