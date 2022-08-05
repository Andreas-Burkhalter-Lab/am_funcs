function fig = ClusterSpikeWfms(wfms,f,t,filters,clustnums,polygons)
% This is the machine-generated representation of a Handle Graphics object
% and its children.  Note that handle values may change when these objects
% are re-created. This may cause problems with any callbacks written to
% depend on the value of the handle at the time the object was saved.
%
% To reopen this object, just type the name of the M-file at the MATLAB
% prompt. The M-file and its associated MAT-file must be on your path.

% Compute their projections on the filters
proj = filters'*wfms;
hfig = Cluster(proj(2,:),proj(1,:));
h1 = uicontrol('Parent',hfig, ...
        'Units','points', ...
        'Position',[412 337 76 30], ...
        'String','AutoCorr',...
        'Callback','ClustSWCallback AutoCorr',...
        'Tag','AutoCorrButton');
h1 = uicontrol('Parent',hfig, ...
        'Units','points', ...
        'Position',[412 298 76 30], ...
        'String','Waveforms',...
        'Callback','ClustSWCallback Waveforms',...
        'Tag','WaveformsButton');
setappdata(hfig,'t',t);
setappdata(hfig,'f',f);
setappdata(hfig,'wfms',wfms);
ymax = max(max(wfms));
ymin = min(min(wfms));
setappdata(hfig,'swylim',[ymin ymax]);
if (nargin > 4)
        % Co-opt the "Clear" function and turn it into a "Revert" function
        h = findobj(hfig,'Tag','ClearButton');
        set(h,'String','Revert','Callback','ClustSPCallback Revert');
        setappdata(hfig,'clustnums0',clustnums);
        setappdata(hfig,'polygons0',polygons);
        ClustSPCallback('Revert',hfig);
end
if nargout > 0, fig = hfig; end
