function result=is_modifier_down(hObject, modifier)
   hFig=get_parent_fig(hObject);
   curModifier=get(hFig, 'selectiontype');
   
   switch lower(modifier)
      case {'s', 'shift'}
	 result=strcmp(curModifier, 'extend');
      case {'c', 'ctrl'}
	 result=strcmp(curModifier, 'alt');
      otherwise
	 error('unsupported keyboard modifier');
   end
