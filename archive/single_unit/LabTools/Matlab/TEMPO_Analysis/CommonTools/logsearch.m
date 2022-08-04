function varargout = LogSearch(varargin)
% LOGSEARCH Application M-file for LogSearch.fig
%    FIG = LOGSEARCH launch LogSearch GUI.
%    LOGSEARCH('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 17-Oct-2001 15:01:25

ProtocolDefs;

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);
    
    % Store a copy of all the keyword arrays, since the sub-functions are retarded.
    handles.dots_keywords = dots_keywords;
    handles.targ_keywords = targ_keywords;
    handles.misc_keywords = misc_keywords;
    handles.neuron_keywords = neuron_keywords;
    handles.one_time_keywords = one_time_keywords;
    handles.eye_calib_keywords = eye_calib_keywords;
    handles.obj_keywords = obj_keywords;
    handles.bar_keywords = bar_keywords;
    handles.bkgnd_keywords = bkgnd_keywords;
    
    % Create 2 arrays, one which holds the internal names of the data sets, and one
    % which holds the user-friendly names of the data sets.
    n = 1;
    for i = 1:NUM_PATCHES
        handles.index_data{n, i} = i;
    end
    handles.internal_dataset{n} = 'handles.dots_keywords'; handles.external_dataset{n} = 'Dots'; n = n + 1;
    for i = 1:NUM_TARGETS
        handles.index_data{n, i} = i;
    end
    handles.internal_dataset{n} = 'handles.targ_keywords'; handles.external_dataset{n} = 'Targets'; n = n + 1;
    handles.internal_dataset{n} = 'handles.misc_keywords'; handles.external_dataset{n} = 'Misc'; n = n + 1;
    handles.internal_dataset{n} = 'handles.neuron_keywords'; handles.external_dataset{n} = 'Neuron'; n = n + 1;
    handles.internal_dataset{n} = 'handles.one_time_keywords'; handles.external_dataset{n} = 'One Time'; n = n + 1;
    handles.internal_dataset{n} = 'handles.eye_calib_keywords'; handles.external_dataset{n} = 'Eye Calib.'; n = n + 1;
    for i = 1:NUM_OBJECTS
        handles.index_data{n, i} = i;
    end
    handles.internal_dataset{n} = 'handles.obj_keywords'; handles.external_dataset{n} = 'Object'; n = n + 1;
    for i = 1:NUM_BARS
        handles.index_data{n, i} = i;
    end
    handles.internal_dataset{n} = 'handles.bar_keywords'; handles.external_dataset{n} = 'Bar'; n = n + 1;
    handles.internal_dataset{n} = 'handles.bkgnd_keywords'; handles.external_dataset{n} = 'Background';  
    
    % Set up a "Truth Matrix", which determine if a variable has already been selected.
    handles.xd = size(handles.internal_dataset, 2);
    handles.yd = size(eval(char(handles.internal_dataset(n))), 2);
    for n = 2:handles.xd
        if handles.yd < size(eval(char(handles.internal_dataset(n))))
            handles.yd = size(eval(char(handles.internal_dataset(n))));
        end
    end
    handles.truthMatrix = [];
    for n = 1:handles.xd
        handles.truthMatrix = [handles.truthMatrix ; 1:handles.yd];
    end
    
    
    % Load the default variables list.
    set(handles.AvailableVarsList, 'String', eval(char(handles.internal_dataset(1))));
    
    % Set up the variable drop down menu and patch drop down.
    set(handles.VarPopUp, 'String', handles.external_dataset);
    set(handles.IndexPopUp, 'String', handles.index_data(1,:));
        
    
    % Save the data.
    guidata(fig, handles);

	if nargout > 0
		varargout{1} = fig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		if (nargout)
			[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
		else
			feval(varargin{:}); % FEVAL switchyard
		end
	catch
		disp(lasterr);
	end

end


%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and 
%| sets objects' callback properties to call them through the FEVAL 
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the 
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.figure1, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.




% --------------------------------------------------------------------
function varargout = DefinedVarsList_Callback(h, eventdata, handles, varargin)
    % Return if nothing valid is selected.
    if isempty(get(h, 'String'));
        return;
    end

    index = get(handles.DefinedVarsList, 'Value');
    answer = handles.storage{index}
    set(handles.RangeEdit1, 'String', answer{1});
    set(handles.RangeEdit2, 'String', answer{2});
    set(handles.NumberEdit, 'String', answer{3});
    set(handles.IndexEdit, 'String', answer{4});
    
    guidata(h, handles);




% --------------------------------------------------------------------
function varargout = VarPopUp_Callback(h, eventdata, handles, varargin)
    % Load the available variables list with the variable set associated with
    % the set picked in the popup menu.
    index = get(h, 'Value');
    set(handles.AvailableVarsList, 'Value', 1);
    set(handles.AvailableVarsList, 'String', eval(char(handles.internal_dataset(index))));
    
    % Set the index PopUp.
    if isempty(handles.index_data(index, 1))
        return;
    end
    
    set(handles.IndexPopUp, 'String', handles.index_data(index,:));
    set(handles.IndexPopUp, 'Value', 1);
    
    guidata(h, handles);




% --------------------------------------------------------------------
function varargout = ChangeVarButton_Callback(h, eventdata, handles, varargin)
    % Return if nothing valid is selected.
    if isempty(get(handles.DefinedVarsList, 'String'));
        return;
    end
    
    answer = inputdlg('Enter new data value: ', 'Change Data Dialog');
    
    % If we get a valid answer from the dialog box, then update the edit box value and
    % the storage data structure.
    if ~isempty(answer)
        set(handles.ValueEdit, 'String', char(answer));
        handles.storage{get(handles.DefinedVarsList, 'Value')} = char(answer);
    end
    
    guidata(h, handles);
    

        
% --------------------------------------------------------------------
function varargout = SetVarButton_Callback(h, eventdata, handles, varargin)
    % Check to see if the item has already been added to the defined list.
    l_index = get(handles.AvailableVarsList, 'Value');
    p_index = get(handles.VarPopUp, 'Value');
    
    
    % Throw up and error if the variable has already been set.
    if handles.truthMatrix(p_index, l_index) == 0
        errordlg('Value already set, Use Change button');
        return;
    end
    
    % Set the flag in th truth matrix.
    handles.truthMatrix(p_index, l_index) = 0;
    
    % Throw up a dialog asking for input.
    answer = inputdlg({'Range Start:', 'Range End:', 'Number:'}, 'Set Variable Dialog');
    
    % Return if the answer was blank.
    if isempty(answer)
        return;
    end
    
    % Put the selected value in the defined variables list.
    l_array = get(handles.AvailableVarsList, 'String');
    set(handles.DefinedVarsList, 'String', [get(handles.DefinedVarsList, 'String') ; l_array(l_index)]);
    set(handles.DefinedVarsList, 'Value', size(get(handles.DefinedVarsList, 'String'), 1));
    
    % Store the location of the variable in the parameters set for easy inclusion into the truth matrix later.
    handles.paramLocation{get(handles.DefinedVarsList, 'Value'), 1, 1} = p_index;
    handles.paramLocation{get(handles.DefinedVarsList, 'Value'), 1, 2} = l_index;
    
    % Store the answer result associated with the variable.
    a = get(handles.DefinedVarsList, 'Value');
    answer{4} = get(handles.IndexPopUp, 'Value');
    handles.storage{a(1)} = answer;
    set(handles.RangeEdit1, 'String', answer{1});
    set(handles.RangeEdit2, 'String', answer{2});
    set(handles.NumberEdit, 'String', answer{3});
    set(handles.IndexEdit, 'String', answer{4});
    
    guidata(h, handles);
    
    


% --------------------------------------------------------------------
function varargout = AvailableVarsList_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = FileMenu_Callback(h, eventdata, handles, varargin)




% ------------------------------------------------------------------------
%   Closes the program when the "File->Exit" item is selected.
%-------------------------------------------------------------------------
function varargout = ExitItem_Callback(h, eventdata, handles, varargin)
    close;


% --------------------------------------------------------------------
function varargout = DelVarButton_Callback(h, eventdata, handles, varargin)
    % Grab String vector and remove the element to be deleted.
    strVector = get(handles.DefinedVarsList, 'String');
    
    % Return if the list is empty.
    if isempty(strVector)
        return;
    end
    
    % Delete the element out of the vector, storage, and the truth matrix.
    index = get(handles.DefinedVarsList, 'Value');
    x = handles.paramLocation{index, 1, 1};
    y = handles.paramLocation{index, 1, 2};
    handles.truthMatrix(x, y) = 1;
    handles.storage(index) = [];
    strVector(index) = [];
    set(handles.DefinedVarsList, 'String', strVector);
    
    % If the vector is empty, then make sure that the edit box shows a blank value.
    if isempty(strVector)
        set(handles.RangeEdit1, 'String', '');
        set(handles.RangeEdit2, 'String', '');
        set(handles.NumberEdit, 'String', '');
        set(handles.IndexEdit, 'String', '');
    else
        set(handles.DefinedVarsList, 'Value', 1);
        answer = handles.storage(1);
        disp(answer{1});
        set(handles.RangeEdit1, 'String', answer{1});
        set(handles.RangeEdit2, 'String', answer{2});
        set(handles.NumberEdit, 'String', answer{3});
        set(handles.IndexEdit, 'String', answer{4});
    end
    
    guidata(h, handles);


% --------------------------------------------------------------------
function varargout = edit4_Callback(h, eventdata, handles, varargin)




% --------------------------------------------------------------------
function varargout = popupmenu3_Callback(h, eventdata, handles, varargin)

