function varargout = cluster_navigator(varargin)
% Usage:
%   cluster_navigator(filename)
%   cluster_navigator(s)
% where
%   filename is the name of the sorting file to open
% OR
%   s is a structure with the following fields:
%     file: the name of the sorting file (i.e., the output of
%       CLUSTER_FITS).
%     progress (default false): if true, shows progress bars

% Copyright 2008 by Timothy E Holy and Zhongsheng Guo
  
%CLUSTER_NAVIGATOR M-file for cluster_navigator.fig
%      CLUSTER_NAVIGATOR, by itself, creates a new CLUSTER_NAVIGATOR or raises the existing
%      singleton*.
%
%      H = CLUSTER_NAVIGATOR returns the handle to a new CLUSTER_NAVIGATOR or the handle to
%      the existing singleton*.
%
%      CLUSTER_NAVIGATOR('Property','Value',...) creates a new CLUSTER_NAVIGATOR using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to cluster_navigator_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      CLUSTER_NAVIGATOR('CALLBACK') and CLUSTER_NAVIGATOR('CALLBACK',hObject,...) call the
%      local function named CALLBACK in CLUSTER_NAVIGATOR.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cluster_navigator

% Last Modified by GUIDE v2.5 30-Jun-2008 16:05:21

% Architectural notes:
%   Each cluster gets a structure, col, that contains handles to
%   already-rendered axes for waveform, autocorrelation, and peak
%   histogram. Panning across the clusters just changes the mapping between
%   the col array and the array of visible axes.
%   Main utility functions:
%     map_get(clusters,clusterID)



% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cluster_navigator_OpeningFcn, ...
                   'gui_OutputFcn',  @cluster_navigator_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

function handlesOut=sortHandleByTags(handles)
   for idx=1:length(handles)
      texpr='^\D+(\d+)$';
      tag=get(handles(idx), 'tag');
      [from, to, tokenExtents, match, tokens]=regexp(tag, texpr, 'once');
      numbers(idx)=str2num(tokens{1});
   end
   [tt, indices]=sort(numbers);
   handlesOut=make_vector(handles(indices), 'row');

   
   
% update [left bottom] only
function updateLB(handle, newLB)
   pos=get(handle, 'position');
   set(handle, 'position', [newLB(1:2) pos(3:4)]);
   
   
function drawWaveform(handles, col)
% Renders the image of the waveform across the array
   if(getappdata(col.waveform, 'valid'))
      return;
   end
   %fprintf('w%s with handle %g\n',col.clusterID,col.waveform);

   setappdata(col.waveform, 'clusterID', col.clusterID);

   hfig=handles.figMain;
   clusters=getappdata(hfig, 'clusters');
   amplitudes=getappdata(hfig, 'amplitudes');
   spiketimes=getappdata(hfig, 'spiketimes');
   residualFile=getappdata(hfig, 'residualFile');
   channels=getappdata(hfig, 'channels');
   sniprange=getappdata(hfig, 'sniprange');
   templatesm=getappdata(hfig, 'templatesm');
   nEventsToPlotAtMost=getappdata(hfig, 'nEventsToPlotAtMost');
   clim=getappdata(hfig, 'waveformClim');
   % memmres = merecmm(replace_parent_dir(residualFile, '')); % TODO: temp
   arrayField = getappdata(hfig,'arrayField');
   
   cluster=map_get(clusters, col.clusterID);
   
   if(isempty(cluster))
      cla(col.waveform); % TODO: even cla(col.waveform, 'reset')?
      return;
   end

   % Get the array geometry
   array_channels = get_hda_holylab(arrayField);

   skip = ceil(length(cluster) / nEventsToPlotAtMost);
   cluster=cluster(1:skip:end);
   % snips = fetch_snippets_from_residual(amplitudes(:,cluster),...
   %    templatesm,spiketimes(cluster),memmres,channels,sniprange);
   snips = templatesm * amplitudes(:,cluster);
   snips = reshape(snips,[size(snips,1)/length(channels) length(channels) size(snips,2)]);
   snips = shape_spikes_to_array(snips,channels,array_channels);

   
   nRows=size(snips, 2);
   nCols=size(snips, 3);
   peakAmp = squeeze(max(abs(mean(snips,4)),[],1));
%    peakAmp=NaN(nRows, nCols);
%    for rowIndex=1:nRows
%       for colIndex=1:nCols
%          tt=squeeze(snips(:, rowIndex, colIndex, :));
%          peakAmp(rowIndex, colIndex)=max(abs(mean(tt, 2)));
%       end
%    end
   
   % axes(col.waveform);
   imagesc(peakAmp, 'parent', col.waveform);
   if(~isempty(clim))
      set(col.waveform, 'clim', clim);
   end
   set(col.waveform,'XTick',[],'YTick',[]);
   
   setappdata(col.waveform, 'valid', 1);
   
   

function drawAutocorr(handles, col)
% Renders the autocorrection histogram
   if(getappdata(col.autocorr, 'valid'))
      return;
   end

   hfig=handles.figMain;
   clusters=getappdata(hfig, 'clusters');
   amplitudes=getappdata(hfig, 'amplitudes');
   spiketimes=getappdata(hfig, 'spiketimes');
   residualFile=getappdata(hfig, 'residualFile');
   channels=getappdata(hfig, 'channels');
   sniprange=getappdata(hfig, 'sniprange');
   templatesm=getappdata(hfig, 'templatesm');
   nEventsToPlotAtMost=getappdata(hfig, 'nEventsToPlotAtMost');
   clim=getappdata(hfig, 'waveformClim');
   scanrate=getappdata(hfig, 'scanrate');
   % memmres = merecmm(replace_parent_dir(residualFile, '')); % TODO: temp
   
   cluster=map_get(clusters, col.clusterID);
   
   if(isempty(cluster))
      cla(col.autocorr);
      return;
   end

   curSpiketimes=spiketimes(cluster);
   trange = getappdata(hfig,'corrTMax');
   tpos=getpixelposition(col.autocorr);
   npix=tpos(3);
   nbins=ceil(npix/2);
   n = autocorrspike(curSpiketimes,trange,nbins);
   binwidth = trange/nbins;
   x = linspace(binwidth/2,trange-binwidth/2,nbins);
   bar(col.autocorr, x,n,'k');
   set(col.autocorr, 'XLim',[0 trange]);
   set(col.autocorr,'XTick',[],'YTick',[]);
   
   setappdata(col.autocorr, 'valid', 1);

   
function drawText(handles, col)
% Draws textual information about each cluster (at the top part of the
% display)
   if(getappdata(col.textclust, 'valid'))
      return;
   end

   hfig=handles.figMain;
   clusters=getappdata(hfig, 'clusters');
   amplitudes=getappdata(hfig, 'amplitudes');
   spiketimes=getappdata(hfig, 'spiketimes');
   residualFile=getappdata(hfig, 'residualFile');
   channels=getappdata(hfig, 'channels');
   sniprange=getappdata(hfig, 'sniprange');
   templatesm=getappdata(hfig, 'templatesm');
   nEventsToPlotAtMost=getappdata(hfig, 'nEventsToPlotAtMost');
   clim=getappdata(hfig, 'waveformClim');
   % memmres = merecmm(replace_parent_dir(residualFile, '')); % TODO: temp
   
   cluster=map_get(clusters, col.clusterID);
   
   if(isempty(cluster))
      n=0;
   else
      n=length(spiketimes(cluster));
   end
   set(col.textclust, 'string', [col.clusterID ':' num2str(n)]);
   
   setappdata(col.textclust, 'valid', 1);
   
   
function drawPeakhist(handles, col)
% Renders the peak-height histogram
   if(getappdata(col.peakhist, 'valid'))
      return;
   end
   
   hfig=handles.figMain;
   clusters=getappdata(hfig, 'clusters');
   amplitudes=getappdata(hfig, 'amplitudes');
   spiketimes=getappdata(hfig, 'spiketimes');
   residualFile=getappdata(hfig, 'residualFile');
   channels=getappdata(hfig, 'channels');
   sniprange=getappdata(hfig, 'sniprange');
   templatesm=getappdata(hfig, 'templatesm');
   nEventsToPlotAtMost=getappdata(hfig, 'nEventsToPlotAtMost');
   clim=getappdata(hfig, 'waveformClim');
   peakVal=getappdata(hfig, 'peakVal');
   thresh=getappdata(hfig, 'thresh');
   
   cluster=map_get(clusters, col.clusterID);
   
   cla(col.peakhist);
   if(isempty(cluster))
      return;
   end
   
   curNormalizedPeakVal=peakVal(cluster)./thresh(cluster)';
   
   % plot(col.peakhist, 5, 5); text(5, 5, ['ph: ' col.clusterID], 'parent', col.peakhist);
   hist(col.peakhist, curNormalizedPeakVal);
   %hold(col.peakhist, 'on');
   ylim=get(col.peakhist, 'ylim');
   hLines=line([-1 1; -1 1], [ylim; ylim]', 'parent', col.peakhist, 'linestyle', '--', 'color', 'r');
   set(col.peakhist,'XTick',[],'YTick',[]);
   
   setappdata(col.peakhist, 'valid', 1);
   

% This function puts the display axes in the appropriate positions on the
% screen, so that the correct ones can be seen.
% re-range the location and redraw if nec.   
function rePosition(handles)
   hfig=handles.figMain;
   onPos=getappdata(hfig, 'onPos');
   offPosLeft=getappdata(hfig, 'offPosLeft');
   offPosRight=getappdata(hfig, 'offPosRight');

   cols=getappdata(hfig, 'cols');
   visibleColStartFrom=getappdata(hfig, 'visibleColStartFrom');
   nVisibleCols=getappdata(hfig, 'nVisibleCols');
   
   offLeft=cols(1:visibleColStartFrom-1);
   onScreen=cols(visibleColStartFrom:min(visibleColStartFrom+nVisibleCols-1, length(cols)));
   offRight=cols(length(offLeft)+length(onScreen)+1:end);
   
   for idx=1:length(offLeft)
      updateLB(offLeft(idx).waveform, offPosLeft);
      updateLB(offLeft(idx).autocorr, offPosLeft);
      updateLB(offLeft(idx).peakhist, offPosLeft);
      updateLB(offLeft(idx).textclust, offPosLeft);
      updateLB(offLeft(idx).checkboxLocked, offPosLeft);
   end
   for idx=1:length(offRight)
      updateLB(offRight(idx).waveform, offPosRight);
      updateLB(offRight(idx).autocorr, offPosRight);
      updateLB(offRight(idx).peakhist, offPosRight);
      updateLB(offRight(idx).textclust, offPosRight);
      updateLB(offRight(idx).checkboxLocked, offPosRight);
   end
   for idx=1:length(onScreen)
      updateLB(onScreen(idx).waveform, onPos.waveform(idx, :));
      updateLB(onScreen(idx).autocorr, onPos.autocorr(idx, :));
      updateLB(onScreen(idx).peakhist, onPos.peakhist(idx, :));
      updateLB(onScreen(idx).textclust, onPos.textclust(idx, :));
      updateLB(onScreen(idx).checkboxLocked, onPos.checkboxLocked(idx, :));
   end
   
   % redraw if nec.
   for idx=1:length(onScreen)
      col=onScreen(idx);
      drawWaveform(handles, col);
      drawAutocorr(handles, col);
      drawText(handles, col);
      drawPeakhist(handles, col);
   end % for, each on-screen col

   
function loadSortResult(hfig, sortResultFile)   
% Loads the fitting results (e.g., amplitudes) for each file
   tt=load(sortResultFile, '-mat');
   spikeClust=tt.spikeClust;
   if(isfield(tt, 'preClusterNavFile'))
      freshStart=false;
      preClusterNavFile=tt.preClusterNavFile;
      clusterIDs=tt.clusterIDs;
      locked = tt.locked;
      tt=load(preClusterNavFile, '-mat');
   else
      freshStart=true;
      preClusterNavFile=sortResultFile;
   end
   
   lminfo=tt.lminfo;
   clust=tt.landmarkClust; % clust(L) is landmark L's clusterID
   spike_labels=spikeClust; % spike_labels(E) is event E's clusterID;
			    % E is also col idx to amplitudes; E is also idx to spiketimes
   if(~freshStart && isequal(clusterIDs{1}, 'noise'))
      spike_labels=spike_labels+1; % due to that agglabel() requires labels>=1
   end% if, has noise
   clusters = agglabel(spike_labels); % clusters{c} is cluster c's events
                                      % that is, clusters{c} is a
                                      % set {E | E is an event for the cluster}
   fitcomp_file=tt.fitFiles;
   n_files = length(fitcomp_file);
   % Load the data
   tstart = zeros(1,n_files);
   nSpikesPerFile = zeros(1,n_files);
   for fileIndex = 1:n_files
     tt=load(fitcomp_file{fileIndex}, '-mat');
     if (fileIndex == 1)
       fine_cluster_svd_file=tt.componentFile;
       % should there be a check that they're all the same??
     end
     h = readheader(tt.fileToFit);
     tstart(fileIndex) = datenum_g([h.date ' ' h.time])*(3600*24);
     amplitudesc{fileIndex}=tt.amplitudes;
     spiketimesc{fileIndex}=tt.spiketimes/h.scanrate;
     nSpikesPerFile(fileIndex) = length(tt.spiketimes);
     peakValc{fileIndex}=tt.peakVal;
     peakChanc{fileIndex}=tt.peakChan;
     merecFile{fileIndex}=tt.fileToFit;
     if isfield(tt,'residualFile')
       residualFile{fileIndex}=tt.residualFile;
     else
       residualFile{fileIndex} = '';
     end
   end
   % Offset the spike times by the file start time
   tstart = tstart-min(tstart);
   for fileIndex = 1:n_files
     spiketimesc{fileIndex} = spiketimesc{fileIndex} + tstart(fileIndex);
     % Do some validation to catch what would otherwise be a hard-to-detect
     % source of crashes in autocorrspike
     if (fileIndex < n_files)
       if (spiketimesc{fileIndex}(end) > tstart(fileIndex+1))
         error('Files appear to be overlapping in time')
       end
     end
   end
   % Concatenate across files
   amplitudes = cat(2,amplitudesc{:});
   n_components = size(amplitudes,1);
   spiketimes = cat(2,spiketimesc{:});
   tstart(end+1) = spiketimes(end);
   peakVal = cat(2,peakValc{:});
   peakChan = cat(2,peakChanc{:});
   
   tt=load(fine_cluster_svd_file, '-mat');
   channels=tt.channels;
   sniprange=tt.sniprange;
   thresh=tt.thresh;
   thresh=thresh(peakChan);
   templates = tt.templates(:,:,1:n_components);
   fu = fit_utilities;
   templatesm = fu.snip2vec_by_t(templates);
   arrayField = tt.arrayField;

   clusterMap=map_new;
   if(freshStart)
      clusterMap=map_put(clusterMap, 'noise', []);
   end
   for clusterIdx=1:length(clusters)
      if(freshStart)
	 clusterID=['c' num2str(clusterIdx)];
      else
	 clusterID=clusterIDs{clusterIdx};
      end
      clusterMap=map_put(clusterMap, clusterID, clusters{clusterIdx});
   end
   

   setappdata(hfig, 'clusters', clusterMap);
   setappdata(hfig, 'amplitudes', amplitudes);
   setappdata(hfig, 'spiketimes', spiketimes);
   setappdata(hfig, 'residualFile', residualFile);
   setappdata(hfig, 'merecFile', merecFile);
   setappdata(hfig, 'channels', channels);
   setappdata(hfig, 'sniprange', sniprange);
   setappdata(hfig, 'templatesm', templatesm);
   setappdata(hfig, 'preClusterNavFile', preClusterNavFile);
   setappdata(hfig, 'peakVal', peakVal);
   setappdata(hfig, 'peakChan', peakChan);   % delete me
   setappdata(hfig, 'thresh', thresh);
   setappdata(hfig, 'arrayField', arrayField);
   setappdata(hfig, 'fileStartTime', tstart);
   setappdata(hfig, 'cumSpikesPerFile', cumsum(nSpikesPerFile));
   if ~freshStart
       setappdata(hfig, 'locked', locked);
   end
   
   
   
function   progressTimerCallback(sender, event_data)
   opt=progress_bar(struct('progress', rand(1)*10, 'max', 10, 'caption', 'Loading ...'));
   figure(opt.handle);

   
% todo: make these func as methods of map
function c=map_new
   c=struct('key', {}, 'value', {});
   
% add a new item or update an existing item   
function c=map_put(c, key, value)
   if(isempty(c))
      c(1).key{1}=key;
      c(1).value{1}=value;
      return;
   end
   t=strmatch(key, c.key, 'exact'); % TODO: use ismemeber() for keys other than strings
   if(isempty(t))
      c.key{end+1}=key;
      c.value{end+1}=value;
   else
      c.key{t}=key;
      c.value{t}=value;
   end
   
function has=map_has(c, key)
   [value, has]=map_get(c, key);
   
function [value, suc]=map_get(c, key)
   value=[];
   suc=0;
   if(isempty(c))
      return;
   end
   t=strmatch(key, c.key, 'exact');
   if(isempty(t))
      return;
   else
      value=c.value{t};
      suc=1;
   end
   
function [c, suc, msg]=map_changekey(c, key, newkey)
   suc=0;
   if(isempty(c))
      error('empty map');
   end
   idx=strmatch(newkey, c.key, 'exact');
   if(~isempty(idx))
      msg='new key exists';
      return;
   end
   idx=strmatch(key, c.key, 'exact');
   if(isempty(idx))
      error('can''t find the key');
   end
   c.key{idx}=newkey;
   
function keys=map_keys(c)
   if(isempty(c))
      keys={};
      return;
   end
   keys=c.key;
   
function [c, suc]=map_erase(c, key)
   suc=0;
   if(isempty(c))
      error('eempty map');
   end
   idx=strmatch(key, c.key, 'exact');
   if(isempty(idx))
      error('can''t find the key');
   end
   c.key(idx)=[];
   c.value(idx)=[];
   
   
function clear_lda_plots(handles)
   hfig=handles.figMain;
   ldaAxes=getappdata(hfig, 'ldaAxes');
   % cla(ldaAxes);      
   for hax=ldaAxes
      cla(hax);
      title(hax,'')
   end
   set(ldaAxes, 'selected', 'off','XTick',[],'YTick',[]);
   
% This updates the "Cluster 1-d overlaps" display when the user
% clicks on a new cluster. We use LDA to get a sense for how different two
% clusters are from each other. 
function update_lda(handles)
   if ~get(handles.checkboxPairwiseUpdate,'Value')
     return
   end
   hfig=handles.figMain;
   clusters=getappdata(hfig, 'clusters');
   amplitudes=getappdata(hfig, 'amplitudes');
   spiketimes=getappdata(hfig, 'spiketimes');
   residualFile=getappdata(hfig, 'residualFile');
   channels=getappdata(hfig, 'channels');
   sniprange=getappdata(hfig, 'sniprange');
   templatesm=getappdata(hfig, 'templatesm');
   
   clear_lda_plots(handles);
   
   selectedWaveform=getappdata(hfig, 'selectedWaveform');
   if(isempty(selectedWaveform))
      return;
   end
   
   max_for_lda = 2000+size(amplitudes,1);  % don't use any more waveforms than this for LDA

   selectedClusterID=getappdata(selectedWaveform, 'clusterID');
   selectedCluster=map_get(clusters, selectedClusterID);
   skip = round(length(selectedCluster)/max_for_lda);
   skip = max(1,skip);
   selectedCluster = selectedCluster(1:skip:end);
   ids=map_keys(clusters);
   usedClusterIDs={};
   tosort=struct('overlap', {[]}, 'proj', {{}}, 'group_label', {{}});
   for idx=1:length(ids)
      clusterID=ids{idx};
      if(ismember(clusterID, {'noise', selectedClusterID}))
         continue;
      end
      cluster=map_get(clusters, clusterID);
      if(isempty(cluster))
         continue;
      end
      usedClusterIDs{end+1}=clusterID;
      
      skip = round(length(cluster)/max_for_lda);
      skip = max(1,skip);
      cluster = cluster(1:skip:end);

      both_events = [selectedCluster, cluster];
      group_label = [ones(size(selectedCluster)), 1+ones(size(cluster))];
      spike_data = amplitudes(:,both_events);
      [eigvec,eigval] = lda(spike_data,group_label);
      proj = eigvec'*spike_data;
      % Quantify the overlap in terms of a d' value
      proj1 = proj(group_label == 1);
      proj2 = proj(group_label == 2);
      dprime = abs(mean(proj1) - mean(proj2)) / sqrt(var(proj1) + var(proj2));
%       n1 = sum(proj1 > min(proj2) & proj1 < max(proj2));
%       n2 = sum(proj2 > min(proj1) & proj2 < max(proj1));
%       tosort.overlap(end+1)=max(n1/length(proj1),n2/length(proj2));
      tosort.overlap(end+1)=dprime;
      tosort.proj{end+1}=proj;
      tosort.group_label{end+1}=group_label;
   end % for, each cluster
   % Put the 
   [sortedEV, indices]=sort(tosort.overlap, 2, 'ascend');
   usedClusterIDs=usedClusterIDs(indices);
   sorted.overlap=tosort.overlap(indices);
   sorted.proj=tosort.proj(indices);
   sorted.group_label=tosort.group_label(indices);
   setappdata(hfig, 'ldaData', sorted);
   setappdata(hfig, 'ldaDisplayFrom', 1);
   setappdata(hfig, 'usedClusterIDs', usedClusterIDs);
   setappdata(hfig, 'selectedLdaID', {});
   
   drawLda(handles);

   
function drawLda(handles)
   hfig=handles.figMain;
   ldaAxes=getappdata(hfig, 'ldaAxes');
   ldaData=getappdata(hfig, 'ldaData');
   ldaDisplayFrom=getappdata(hfig, 'ldaDisplayFrom');
   usedClusterIDs=getappdata(hfig, 'usedClusterIDs');
   nLdaToDraw=min(length(usedClusterIDs)-ldaDisplayFrom+1, length(ldaAxes));
   logScaled = get(handles.checkboxLog,'Value');
   
   for idx=1:nLdaToDraw
      hax=ldaAxes(idx);
      % axes(hax);
      cla(hax);
      thisPairIndex = ldaDisplayFrom+idx-1;
      proj=ldaData.proj{thisPairIndex};
      group_label=ldaData.group_label{thisPairIndex};
      clusterID=usedClusterIDs{thisPairIndex};
      overlap=ldaData.overlap(thisPairIndex);
      
      %set(hax,'YScale','linear');
      hold(hax, 'on');
      [y,x] = hist(proj(group_label == 2),20);
      %hb = bar(x,y+1,'b');
      %set(hb,'EdgeColor','none');
      plot(hax, x,y+1,'b');
      [y,x] = hist(proj(group_label == 1),20);
      %hb = bar(x,y+1,'r');
      %set(hb,'EdgeColor','none');
      plot(hax, x,y+1,'r');
      if logScaled
        set(hax,'YScale','log')
      else
        set(hax,'YScale','linear')
      end
      set(hax,'XTick',[],'YTickMode','auto');
      title(hax, [clusterID ': ' sprintf('%0.2g',overlap)]);
   end

function handled=ldaAxes_mouse_up(sender, event_args)
   handled = 1;
   hax=sender;
   handles=guidata(hax);
   hfig=handles.figMain;
   selectedLdaID=getappdata(hfig, 'selectedLdaID');
   
   if(isequal('on', get(hax, 'selected')))
      set(hax, 'selected', 'off');
      str2del=get(get(hax, 'title'), 'string');
      selectedLdaID=setdiff(selectedLdaID, {str2del});
   else
      clusterID=get(get(hax, 'title'), 'string');
      cIndex = strfind(clusterID,':');
      clusterID = clusterID(1:cIndex-1);
      cols=getappdata(hfig, 'cols');
      idx=strmatch(clusterID, {cols.clusterID}, 'exact');
      col=cols(idx);
      hcheckbox=col.checkboxLocked;
      checked=get(hcheckbox,'Value') == get(hcheckbox,'Max');
      if(checked)
	 btn=questdlg('The cluster is locked, do you want to select it?', ...
		      'Cluster is locked', 'yes', 'no', 'no');
	 if(isequal(btn, 'no'))
	    return;
	 end
      end
      
      set(hax, 'selected', 'on');
      selectedLdaID{end+1}=clusterID;
   end
   setappdata(hfig, 'selectedLdaID', selectedLdaID);
   
   
function handled=waveformAxes_mouse_up(sender, event_args)
   handled = 1;
   hax=sender;
   handles=guidata(hax);
   hfig=handles.figMain;
   oldSelection=getappdata(hfig, 'selectedWaveform');
   
   if(oldSelection==hax)
      return;
   end
   
   clusters=getappdata(hfig, 'clusters');
   specialClusters=getappdata(hfig, 'specialClusters');
   amplitudes=getappdata(hfig, 'amplitudes');
   spiketimes=getappdata(hfig, 'spiketimes');
   residualFile=getappdata(hfig, 'residualFile');
   channels=getappdata(hfig, 'channels');
   sniprange=getappdata(hfig, 'sniprange');
   templatesm=getappdata(hfig, 'templatesm');

   id=getappdata(hax, 'clusterID');
   if(isequal(id, 'noise'))
      return; 
   end
   if(isempty(map_get(clusters, id)))
      return;
   end
   
   set(oldSelection, 'selected', 'off');
   set(hax, 'selected', 'on');
   setappdata(hfig, 'selectedWaveform', hax);

   update_lda(handles);

   
   
function invalidateCol(handles, clusterID)
   hfig=handles.figMain;
   cols=getappdata(hfig, 'cols');
   idx=strmatch(clusterID, {cols.clusterID}, 'exact');
   col=cols(idx);
   for h=[col.waveform col.autocorr col.peakhist col.textclust col.checkboxLocked]
      setappdata(h, 'valid', 0);
   end
   
   
function stickCol(sender, event_args)
   hfig=sender;
   handles=guidata(hfig);
   selectedWaveform=getappdata(hfig, 'selectedWaveform');
   if(isempty(selectedWaveform))
      return; % todo: may need warn user about this
   end
   clusterID=getappdata(selectedWaveform, 'clusterID');
   cols=getappdata(hfig, 'cols');
   idx=strmatch(clusterID, {cols.clusterID}, 'exact');
   visibleColStartFrom=getappdata(hfig, 'visibleColStartFrom');
   nVisibleCols=getappdata(hfig, 'nVisibleCols');
   if(idx<visibleColStartFrom || idx>=visibleColStartFrom+nVisibleCols)
      return; % todo: may need warn user
   end
   
   sticky=getappdata(hfig, 'sticky');
   sticky(idx-visibleColStartFrom+1) = ~sticky(idx-visibleColStartFrom+1);
   setappdata(hfig, 'sticky', sticky);
   
   
   
function delCol(sender, event_args)   
   hfig=sender;
   handles=guidata(hfig);
   selectedWaveform=getappdata(hfig, 'selectedWaveform');
   if(isempty(selectedWaveform))
      return; % todo: may need warn user about this
   end
   
   % clear the selection
   setappdata(hfig, 'selectedWaveform', []);
   set(selectedWaveform, 'selected', 'off');
   update_lda(handles);

   selectedClusterID=getappdata(selectedWaveform, 'clusterID');
   cols=getappdata(hfig, 'cols');
   clusters=getappdata(hfig, 'clusters');
   idx=find([cols.waveform]==selectedWaveform);
   col=cols(idx);
   cols(idx)=[];
   uninstall_mouse_event_handler(col.waveform, 'up', @waveformAxes_mouse_up);
   arrayfun(@delete, [col.waveform col.autocorr col.peakhist col.textclust col.checkboxLocked]);
   clusterToDel=map_get(clusters, selectedClusterID);
   clusters=map_erase(clusters, selectedClusterID);
   noise=map_get(clusters, 'noise');
   noise=sort([noise clusterToDel]);
   clusters=map_put(clusters, 'noise', noise);
   
   setappdata(hfig, 'cols', cols);
   setappdata(hfig, 'clusters', clusters);
   
   invalidateCol(handles, 'noise');
   
   rePosition(handles);
   
   
   
% --- Executes just before cluster_navigator is made visible.
function cluster_navigator_OpeningFcn(hObject, eventdata, handles, varargin)
   % This function has no output args, see OutputFcn.
   % hObject    handle to figure
   % eventdata  reserved - to be defined in a future version of MATLAB
   % handles    structure with handles and user data (see GUIDATA)
   % varargin   unrecognized PropertyName/PropertyValue pairs from the
   %            command line (see VARARGIN)
   
   % Choose default command line output for cluster_navigator
   handles.output = hObject;
   
   if isstruct(varargin{1})
     navParams=varargin{1};
   else
     navParams.file = varargin{1};
   end
   navParams = default(navParams,'progress',false);
   
   if navParams.progress
     progtimer=timer('TimerFcn', @progressTimerCallback, ...
       'ExecutionMode', 'fixedRate', 'period', 0.1);
     start(progtimer);
   end
   
   % figure out placeholders's position
   hfig=hObject;
   waveformAxes=sortHandleByTags(findobj(hfig, '-regexp', 'tag', 'axesW.+'));
   autocorrAxes=sortHandleByTags(findobj(hfig, '-regexp', 'tag', 'axesAC.+'));
   peakhistAxes=sortHandleByTags(findobj(hfig, '-regexp', 'tag', 'axesPh.+'));
   textclusts=sortHandleByTags(findobj(hfig, '-regexp', 'tag', 'textClust.+'));
   checkboxLockeds=sortHandleByTags(findobj(hfig, '-regexp', 'tag', 'checkboxLocked.+'));
   
   ldaAxes=sortHandleByTags(findobj(hfig, '-regexp', 'tag', 'axesOl.+'));
   
   % TODO: temp
   % set(get(hfig, 'children'), 'units', 'norm');
   
   % on screen positions
   onPos.waveform=cell2mat(get(waveformAxes, 'position'));
   onPos.autocorr=cell2mat(get(autocorrAxes, 'position'));
   onPos.peakhist=cell2mat(get(peakhistAxes, 'position'));
   onPos.textclust=cell2mat(get(textclusts, 'position'));
   onPos.checkboxLocked=cell2mat(get(checkboxLockeds, 'position'));
   
   % off screen pos 
   offPosLeft=[-1 -1]*3000;
   offPosRight=[2 2]*3000;
   
   sortResultFile=navParams.file;
   loadSortResult(hfig, sortResultFile);
   clusters=getappdata(hfig, 'clusters');
   clusterIDs=map_keys(clusters);
   nClusters=length(clusterIDs);
   
   
   offLeft=struct('clusterID', {}, 'waveform', {}, 'autocorr', {}, ...
      'peakhist', {}, 'textclust', {}, 'checkboxLocked', {});
   nVisibleCols=length(waveformAxes);
   nOn=min(nVisibleCols, nClusters);
   nOffRight=max(nClusters-nVisibleCols, 0);
   
   sticky=false(1, nOn);
   %sticky([1 5])=true; % TODO: make it just 1
   
   for idx=1:nOn
      onScreen(idx).clusterID=clusterIDs{idx};
      onScreen(idx).waveform=copyobj(waveformAxes(idx), hfig);
      plot(onScreen(idx).waveform, 5, 5); text(5, 5, ['w ' num2str(idx)], 'parent', onScreen(idx).waveform);
      onScreen(idx).autocorr=copyobj(autocorrAxes(idx), hfig);
      plot(onScreen(idx).autocorr, 5, 5); text(5, 5, ['ac ' num2str(idx)], 'parent', onScreen(idx).autocorr);
      onScreen(idx).peakhist=copyobj(peakhistAxes(idx), hfig);
      plot(onScreen(idx).peakhist, 5, 5); text(5, 5, ['ph ' num2str(idx)], 'parent', onScreen(idx).peakhist);
      onScreen(idx).textclust=copyobj(textclusts(idx), hfig);
      set(onScreen(idx).textclust, 'string', [num2str(idx) ': 0']); 
      onScreen(idx).checkboxLocked=copyobj(checkboxLockeds(idx), hfig);
      % set(onScreen(idx).checkboxLocked, 'tooltipstring', num2str(idx));
      
      setappdata(onScreen(idx).waveform, 'valid', 0); % i.e. invalidated, need redrawing
      % setappdata(onScreen(idx).waveform, 'clusterID', onScreen(idx).clusterID);
   end

   offRight=offLeft;
   for idx=1:nOffRight
      idxAll=idx+nOn;
      offRight(idx).clusterID=clusterIDs{idxAll}; 
      offRight(idx).waveform=copyobj(waveformAxes(1), hfig);
      plot(offRight(idx).waveform, 5, 5); text(5, 5, ['w ' num2str(idxAll)], 'parent', offRight(idx).waveform);
      offRight(idx).autocorr=copyobj(autocorrAxes(1), hfig);
      plot(offRight(idx).autocorr, 5, 5); text(5, 5, ['ac ' num2str(idxAll)], 'parent', offRight(idx).autocorr);
      offRight(idx).peakhist=copyobj(peakhistAxes(1), hfig);
      plot(offRight(idx).peakhist, 5, 5); text(5, 5, ['ph ' num2str(idxAll)], 'parent', offRight(idx).peakhist);
      offRight(idx).textclust=copyobj(textclusts(1), hfig);
      set(offRight(idx).textclust, 'string', [num2str(idxAll) ': 0']); 
      offRight(idx).checkboxLocked=copyobj(checkboxLockeds(1), hfig);
      % set(offRight(idx).checkboxLocked, 'tooltipstring', num2str(idxAll));
   
      setappdata(offRight(idx).waveform, 'valid', 0); % i.e. invalidated, need redrawing
      % setappdata(offRight(idx).waveform, 'clusterID', offRight(idx).clusterID);
   end
   
   cols=[offLeft onScreen offRight];
   if isappdata(hfig,'locked')
       locked = getappdata(hfig,'locked');
       hcb = [cols.checkboxLocked];
       set(hcb,{'Value'},mat2cell(locked(:),ones(1,length(locked)),1));
   end
   
   set([cols.waveform cols.autocorr cols.peakhist], ...
      'ylim', [0 nClusters+1], 'ytick', [], 'xtick', []);
   
   set(hfig, 'units', 'pixels'); % TODO: tmp
   allWaveformAxes=[cols.waveform];
   for hax=allWaveformAxes
      install_mouse_event_handler(hax, 'up', @waveformAxes_mouse_up);
   end
   for hax=ldaAxes
      install_mouse_event_handler(hax, 'up', @ldaAxes_mouse_up);
   end
   
   bind_shortcut(hfig, 'd', @delCol);
   bind_shortcut(hfig, 'backspace', @delCol);
   bind_shortcut(hfig, 'delete', @delCol);
   bind_shortcut(hfig, 's', @stickCol);
   
   placeHolders={waveformAxes, autocorrAxes, peakhistAxes, textclusts, checkboxLockeds};
   set(cat(2, placeHolders{:}), ...
      'visible', 'off');
   
   setappdata(hfig, 'onPos', onPos);
   setappdata(hfig, 'offPosLeft', offPosLeft);
   setappdata(hfig, 'offPosRight', offPosRight);
   setappdata(hfig, 'cols', cols);
   setappdata(hfig, 'visibleColStartFrom', 1);
   setappdata(hfig, 'nVisibleCols', nVisibleCols);
   setappdata(hfig, 'sticky', sticky);
   setappdata(hfig, 'nEventsToPlotAtMost', 60);
   setappdata(hfig, 'ldaAxes', ldaAxes);
   setappdata(hfig, 'corrTMax', 0.5);
   
   % Update handles structure
   guidata(hObject, handles);
   
   % precomp so we can a global clim
   for idx=1:length(cols)
      col=cols(idx);
      drawWaveform(handles, col);
   end % for, each on-screen col
   haxes=[cols.waveform];
   clims=get(haxes, 'clim');
   clims=cell2mat(clims);
   clim=[0 max(clims(:,2))];
   set(haxes, 'clim', clim);
   setappdata(hfig, 'waveformClim', clim);
   
   rePosition(handles);

   if navParams.progress
     stop(progtimer);
     delete(progtimer);
     progress_bar(struct('progress', 10, 'max', 10));
   end
   
   % UIWAIT makes cluster_navigator wait for user response (see UIRESUME)
   % uiwait(handles.figMain);
   

% --- Outputs from this function are returned to the command line.
function varargout = cluster_navigator_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function shiftColOnce(handles, shiftToRight)
   hfig=handles.figMain;
   cols=getappdata(hfig, 'cols');
   visibleColStartFrom=getappdata(hfig, 'visibleColStartFrom');
   nVisibleCols=getappdata(hfig, 'nVisibleCols');
   sticky=getappdata(hfig, 'sticky');
   
   offLeft=cols(1:visibleColStartFrom-1);
   onScreen=cols(visibleColStartFrom:min(visibleColStartFrom+nVisibleCols-1, length(cols)));
   offRight=cols(length(offLeft)+length(onScreen)+1:end);
   sticky = sticky(1:length(onScreen));
   
   if(shiftToRight)
      % shift to right
      if(isempty(offLeft)) return; end
      candidates=[offLeft(end) onScreen(~sticky)];
      offRight=[candidates(end) offRight];
      onScreen(~sticky)=candidates(1:end-1);
      offLeft=offLeft(1:end-1);
   else
      % shift to left
      if(isempty(offRight)) return; end
      candidates=[onScreen(~sticky) offRight(1)];
      offLeft(end+1)=candidates(1);
      onScreen(~sticky)=candidates(2:end);
      offRight=offRight(2:end);
   end
   
   cols = [offLeft onScreen offRight];
   setappdata(hfig, 'cols', cols);
   setappdata(hfig, 'visibleColStartFrom', length(offLeft)+1);

   
function shiftCol(handles, nCol)
   shiftToRight=nCol>0;
   nCol=abs(nCol);
   for i=1:nCol
      shiftColOnce(handles, shiftToRight);
   end
   rePosition(handles);
      
% --- Executes on button press in btnShow1moreLeft.
function btnShow1moreLeft_Callback(hObject, eventdata, handles)
% hObject    handle to btnShow1moreLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   shiftCol(handles, 1);

% --- Executes on button press in btnShow1moreRight.
function btnShow1moreRight_Callback(hObject, eventdata, handles)
% hObject    handle to btnShow1moreRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   shiftCol(handles, -1);

% --- Executes on button press in btnShow5moreRight.
function btnShow5moreRight_Callback(hObject, eventdata, handles)
% hObject    handle to btnShow5moreRight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   shiftCol(handles, -5);

% --- Executes on button press in btnShow5moreLeft.
function btnShow5moreLeft_Callback(hObject, eventdata, handles)
% hObject    handle to btnShow5moreLeft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   shiftCol(handles, 5);

% --- Executes on button press in btnShowLeftmost.
function btnShowLeftmost_Callback(hObject, eventdata, handles)
% hObject    handle to btnShowLeftmost (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   hfig=handles.figMain;
   visibleColStartFrom=getappdata(hfig, 'visibleColStartFrom');
   shiftCol(handles, visibleColStartFrom-1);   

% --- Executes on button press in btnShowRightmost.
function btnShowRightmost_Callback(hObject, eventdata, handles)
% hObject    handle to btnShowRightmost (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   hfig=handles.figMain;

   cols=getappdata(hfig, 'cols');
   visibleColStartFrom=getappdata(hfig, 'visibleColStartFrom');
   nVisibleCols=getappdata(hfig, 'nVisibleCols');
   
   offLeft=cols(1:visibleColStartFrom-1);
   onScreen=cols(visibleColStartFrom:min(visibleColStartFrom+nVisibleCols-1, length(cols)));
   offRight=cols(length(offLeft)+length(onScreen)+1:end);
   
   shiftCol(handles, -length(offRight));   


function shiftLda(handles, nLda)
   hfig=handles.figMain;
   ldaAxes=getappdata(hfig, 'ldaAxes');
   ldaDisplayFrom=getappdata(hfig, 'ldaDisplayFrom');
   usedClusterIDs=getappdata(hfig, 'usedClusterIDs');
   oldLdaDisplayFrom=ldaDisplayFrom;
   ldaDisplayFrom=ldaDisplayFrom-nLda;
   if(ldaDisplayFrom<0)
      ldaDisplayFrom=1;
   end
   maxFrom=max(1, length(usedClusterIDs)-length(ldaAxes)+1);
   if(ldaDisplayFrom>maxFrom)
      ldaDisplayFrom=maxFrom;
   end
   
   if(ldaDisplayFrom~=oldLdaDisplayFrom)
      setappdata(hfig, 'ldaDisplayFrom', ldaDisplayFrom);
      drawLda(handles);
   end
   
% --- Executes on button press in btnShow1moreLeftLda.
function btnShow1moreLeftLda_Callback(hObject, eventdata, handles)
% hObject    handle to btnShow1moreLeftLda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   shiftLda(handles, 1);

% --- Executes on button press in btnShow1moreRightLda.
function btnShow1moreRightLda_Callback(hObject, eventdata, handles)
% hObject    handle to btnShow1moreRightLda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   shiftLda(handles, -1);

% --- Executes on button press in btnShow5moreRightLda.
function btnShow5moreRightLda_Callback(hObject, eventdata, handles)
% hObject    handle to btnShow5moreRightLda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   shiftLda(handles, -5);

% --- Executes on button press in btnShow5moreLeftLda.
function btnShow5moreLeftLda_Callback(hObject, eventdata, handles)
% hObject    handle to btnShow5moreLeftLda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   shiftLda(handles, 5);

% --- Executes on button press in btnShowLeftmostLda.
function btnShowLeftmostLda_Callback(hObject, eventdata, handles)
% hObject    handle to btnShowLeftmostLda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   shiftLda(handles, 5000);

% --- Executes on button press in btnShowRightmostLda.
function btnShowRightmostLda_Callback(hObject, eventdata, handles)
% hObject    handle to btnShowRightmostLda (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   shiftLda(handles, -5000);

% --- Executes on button press in checkboxLocked2.
function checkboxLocked2_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxLocked2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxLocked2


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on button press in checkboxLocked6.
function checkboxLocked6_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxLocked6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxLocked6


% --- Executes on button press in checkboxLocked5.
function checkboxLocked5_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxLocked5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxLocked5


% --- Executes on button press in checkboxLocked4.
function checkboxLocked4_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxLocked4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxLocked4


% --- Executes on button press in checkboxLocked7.
function checkboxLocked7_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxLocked7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxLocked7


% --- Executes on button press in checkboxLocked10.
function checkboxLocked10_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxLocked10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxLocked10


% --- Executes on button press in checkboxLocked9.
function checkboxLocked9_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxLocked9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxLocked9


% --- Executes on button press in checkboxLocked8.
function checkboxLocked8_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxLocked8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxLocked8

function ids=createNewClusterIDs(handles, n)
   if (n < 1)
     ids = {};
     return;
   end
   hfig=handles.figMain;
   clusters=getappdata(hfig, 'clusters');
   clusterIDs=map_keys(clusters);
   nums=zeros(1, length(clusterIDs));
   for idx=1:length(nums)
      t=str2num(clusterIDs{idx}(2:end));
      if(~isempty(t))
	 nums(idx)=t(1);
      end
   end
   maxnum=max([nums 1000]);
   newnums=maxnum+(1:n);
   ids=cellfun(@(x) (['c' num2str(x)]), num2cell(newnums), 'UniformOutput', false);
   
   

% --- Executes on button press in btnRecluster.
function btnRecluster_Callback(hObject, eventdata, handles)
% hObject    handle to btnRecluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   hfig=handles.figMain;
   selectedLdaID=getappdata(hfig, 'selectedLdaID');
   
   selectedWaveform=getappdata(hfig, 'selectedWaveform');
   if(isempty(selectedWaveform))
      return;
   end
   selectedClusterID=getappdata(selectedWaveform, 'clusterID');
   reClusterIDs=[{selectedClusterID} selectedLdaID];
   
   clusters=getappdata(hfig, 'clusters');
   amplitudes=getappdata(hfig, 'amplitudes');
   spiketimes=getappdata(hfig, 'spiketimes');
   tstart = getappdata(hfig,'fileStartTime');
   corrTMax = getappdata(hfig,'corrTMax');
   bigcluster=[];
   labels=[];
   for idx=1:length(reClusterIDs)
      cluster=map_get(clusters, reClusterIDs{idx});
      bigcluster=[bigcluster cluster];
      labels    =[labels idx*ones(1, length(cluster))];
   end
   points=amplitudes(:, bigcluster);
   mdcluster_dataops = struct('points', points, 'time', spiketimes(bigcluster), 'timerange', tstart([1 end]), 'corrTMax', corrTMax);
   if (length(reClusterIDs) > 1)
     mdcluster_dataops.labels = labels;
   end
   mdcluster_methodops = mdcluster_options(struct('use_maxpoints_in_projection',true,'maxpoints',4000,'default_n_dims',3));
   if (length(bigcluster) > 5000)
     mdcluster_methodops.markerSizes = [4 6];
   end
   if (length(reClusterIDs) == 1)
     mdcluster_methodops.defaultProjectionMethod = 'PCA';
   else
     mdcluster_methodops.defaultProjectionMethod = 'LDA';
     ldapmIndex = strmatch('LDA',{mdcluster_methodops.projectionMethods.name},'exact');
     tmp = mdcluster_methodops.projectionMethods(ldapmIndex);
     ldapmPIndex = strmatch('# dimensions',{tmp.parameterList.name},'exact');
     if ~isempty(ldapmPIndex)
       tmp.parameterList(ldapmPIndex).value = length(reClusterIDs)-1;
     end
     mdcluster_methodops.projectionMethods(ldapmIndex) = tmp;
   end
   
   % Do the new clustering
   newlabels=mdcluster(mdcluster_dataops,mdcluster_methodops);
   isCancelled=isempty(newlabels);
   if(isCancelled)
      return;
   end
   
   % clear the selection
   setappdata(hfig, 'selectedWaveform', []);
   set(selectedWaveform, 'selected', 'off');
   update_lda(handles);
   
   % Convert the new cluster numbers into an index into the collection of
   % all spikes
   [newClusters,nSpikesPerNewCluster]=agglabel(newlabels);
   nNewClusters=length(newClusters);
   for label=1:nNewClusters
      newClusters{label}=sort(bigcluster(newClusters{label}));
   end
   % Sort them so that the clusters with the fewest spikes appear earliest
   % (that way we get consistent behavior no matter what cluster #s are
   % assigned)
   [tmp,sortIndex] = sort(nSpikesPerNewCluster);
   newClusters = newClusters(sortIndex);
   
   % Update the information for user display; we'll use new IDs even if we
   % are replacing, so that the "identity" is marked as being different
   cols=getappdata(hfig, 'cols');
   onPos=getappdata(hfig, 'onPos');
   newClusterIDs=createNewClusterIDs(handles, length(newClusters));
   nReplaced = length(reClusterIDs);
   nReplaced = 0;
   nNew = length(newClusters)-nReplaced;
   delIndex = 1:length(reClusterIDs);

   % Do the replacements
   for idx=1:min(nReplaced,length(newClusters))
      clusterID=reClusterIDs{idx};
      colIdx=strmatch(clusterID, {cols.clusterID}, 'exact');
      cols(colIdx).clusterID=newClusterIDs{idx};
      % updatedColIndices(end+1)=colIdx;
      clusters=map_put(clusters, clusterID, newClusters{idx});
      clusters=map_changekey(clusters, clusterID, cols(colIdx).clusterID); 
   end
   

   % Add new cols
   newcols=struct('clusterID', {}, 'waveform', {}, 'autocorr', {}, ...
      'peakhist', {}, 'textclust', {}, 'checkboxLocked', {});
   for idx=1:nNew
      newcols(idx).clusterID=newClusterIDs{idx+nReplaced};
      newcols(idx).waveform=axes('Parent',hfig,'Units','pixels','position', onPos.waveform(1,:));
      newcols(idx).autocorr=axes('Parent',hfig,'Units','pixels','position', onPos.autocorr(1,:));
      newcols(idx).peakhist=axes('Parent',hfig,'Units','pixels','position', onPos.peakhist(1,:));
      newcols(idx).textclust=copyobj(cols(1).textclust, hfig);
      newcols(idx).checkboxLocked=copyobj(cols(1).checkboxLocked, hfig);
      install_mouse_event_handler(newcols(idx).waveform, 'up', @waveformAxes_mouse_up);
      
      clusters=map_put(clusters, newcols(idx).clusterID, ...
		       newClusters{idx+nReplaced});
   end
   
   % deletion of some cols
   for idx=delIndex
      %clusterID=reClusterIDs{length(newClusters)+idx};
      clusterID=reClusterIDs{idx};
      colIdx=strmatch(clusterID, {cols.clusterID}, 'exact');
      col=cols(colIdx);
      cols(colIdx)=[];
      uninstall_mouse_event_handler(col.waveform, 'up', @waveformAxes_mouse_up);
      arrayfun(@delete, [col.waveform col.autocorr col.peakhist col.textclust col.checkboxLocked]);
      clusters=map_erase(clusters, clusterID);
   end

   setappdata(hfig, 'cols', [cols newcols]);
   setappdata(hfig, 'clusters', clusters);
   
   for idx=1:length(newClusterIDs)
      invalidateCol(handles, newClusterIDs{idx});
   end
   
   rePosition(handles);

   
   
% --- Executes on button press in buttonTimeMarkers.
function buttonTimeMarkers_Callback(hObject, eventdata, handles)
% hObject    handle to buttonTimeMarkers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btnViewWaveforms.
function btnViewWaveforms_Callback(hObject, eventdata, handles)
% hObject    handle to btnViewWaveforms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   hfig=handles.figMain;
   selectedWaveform=getappdata(hfig, 'selectedWaveform');
   if(isempty(selectedWaveform))
      return; % todo: may need warn user about this
   end
   selectedClusterID=getappdata(selectedWaveform, 'clusterID');
   clusters=getappdata(hfig, 'clusters');
   amplitudes=getappdata(hfig, 'amplitudes');
   spiketimes=getappdata(hfig, 'spiketimes');
   residualFile=getappdata(hfig, 'residualFile');
   channels=getappdata(hfig, 'channels');
   sniprange=getappdata(hfig, 'sniprange');
   templatesm=getappdata(hfig, 'templatesm');
   nEventsToPlotAtMost=getappdata(hfig, 'nEventsToPlotAtMost');
   clim=getappdata(hfig, 'waveformClim');

   cluster=map_get(clusters, selectedClusterID);
   if(isempty(cluster))
      return; % TODO: may want to warn user
   end
   skip = ceil(length(cluster) / nEventsToPlotAtMost);
   if (skip < 1)
     skip = 1;
   end
   cluster=cluster(1:skip:end);
   snips = templatesm * amplitudes(:,cluster);
   snips = reshape(snips,[size(snips,1)/length(channels) length(channels) size(snips,2)]);
   arrayField = getappdata(hfig,'arrayField');
   array_channels = get_hda_holylab(arrayField);
   snips = shape_spikes_to_array(snips,channels,array_channels);
   figure;
   hax = SplitVert([0.2 0.25],[1 0 1]);
   haxh = SplitHoriz([0.6 0.7],[1 0 1],hax(2));
   plot_spikes_on_array(snips,hax(1),struct('showscale',true));
   plot(haxh(1),amplitudes(:,cluster));
   peakChan = getappdata(hfig,'peakChan');
   peakVal = getappdata(hfig,'peakVal');
   plot(haxh(2),peakChan(cluster),peakVal(cluster),'.');
   uicontrol('Style','pushbutton','String','Fetch',...
     'Position',[0 0 50 25],'Callback',{@clustnav_fetch,handles,cluster});
   
function clustnav_fetch(hObject,eventData,handles,cluster)
hfig = handles.figMain;
spiketimes = getappdata(hfig, 'spiketimes');
cumSpikesPerFile = [0 getappdata(hfig, 'cumSpikesPerFile')];
merecFile = getappdata(hfig, 'merecFile');
channels = getappdata(hfig, 'channels');
sniprange = getappdata(hfig, 'sniprange');
fileStartTime = getappdata(hfig, 'fileStartTime');
arrayField = getappdata(hfig,'arrayField');
peakChan = getappdata(hfig,'peakChan');
peakVal = getappdata(hfig,'peakVal');
% Parcel the events back in to files, and load the snippet
nFiles = length(merecFile);
spikeTimeFile = cell(1,nFiles);
snips = {};
for i = 1:nFiles
  c = cluster(cluster > cumSpikesPerFile(i) & cluster < cumSpikesPerFile(i+1));
  t = spiketimes(c) - fileStartTime(i);
  memm = merecmm(merecFile{i});
  tscan = round(t*memm.scanrate);
  for j = 1:length(tscan)
    snips{end+1} = memm(channels,tscan(j)+(sniprange(1):sniprange(2)))';
  end
end
peakChan = peakChan(cluster);
peakVal = peakVal(cluster);
% Plot it
snips = cat(3,snips{:});
array_channels = get_hda_holylab(arrayField);
snips = shape_spikes_to_array(snips,channels,array_channels);
pVArray = zeros(length(channels),length(cluster));
indx = sub2ind(size(pVArray),peakChan,1:length(cluster));
pVArray(indx) = peakVal;
pVArray = shape_spikes_to_array(shiftdim(pVArray,-1),channels,array_channels);
figure
plot_spikes_on_array(snips,struct('showscale',true));
% hax = get(gcf,'Children');
% set(hax,'NextPlot','add')
% figure
% hline = plot_spikes_on_array(pVArray,struct('showscale',true));
% set(hline,'Marker','x')


% --- Executes on button press in btnSave.
function btnSave_Callback(hObject, eventdata, handles)
% hObject    handle to btnSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   hfig=handles.figMain;
   clusterNavFile=getappdata(hfig, 'clusterNavFile');
   if(isempty(clusterNavFile))
      btnSaveAs_Callback(hObject, eventdata, handles);
      return;
   end
   preClusterNavFile=getappdata(hfig, 'preClusterNavFile');
   clusters=getappdata(hfig, 'clusters');
   ids=map_keys(clusters);
   intLabel=0;
   bigcluster=[]; % event times
   spikeClust=[]; % the intergral labelling
   clusterIDs={}; % reordered clusterIDs (noise cames at the first)
   t=strmatch('noise', ids, 'exact'); % NOTE: 'noise' clusters may be renamed
   if(~isempty(t))
      % move noise to the beginning of ids
      ids(t)=[];
      ids=[{'noise'} ids];
   else
      intLabel=1; % no noise cluster
   end
   
   for clusterIdx=1:length(ids)
      clusterID=ids{clusterIdx};
      cluster=map_get(clusters, clusterID);
      bigcluster=[bigcluster cluster];
      spikeClust=[spikeClust intLabel*ones(1, length(cluster))];
      clusterIDs{end+1}=clusterID;
      intLabel=intLabel+1;
   end
   cols=getappdata(hfig, 'cols');
   lockedh = [cols.checkboxLocked];
   locked = cell2mat(get(lockedh,'Value'))';
   locked = locked(findainb(clusterIDs,{cols.clusterID}));
   [tt, indices]=sort(bigcluster);
   spikeClust=spikeClust(indices);
   save(clusterNavFile, 'preClusterNavFile', ...
	'spikeClust', ...
	'clusterIDs', ...
  'locked' ...
	);
   

% --- Executes on button press in btnQuit.
function btnQuit_Callback(hObject, eventdata, handles)
% hObject    handle to btnQuit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   close(handles.figMain);

% --- Executes on button press in checkboxLocked1.
function checkboxLocked1_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxLocked1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxLocked1


% --- Executes on button press in btnSaveAs.
function btnSaveAs_Callback(hObject, eventdata, handles)
% hObject    handle to btnSaveAs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
   hfig=handles.figMain;
   [filename, pathname] = uiputfile('*.cluster_nav', 'Pick a cluster navigator result file');
   if isequal(filename,0) || isequal(pathname,0)
      return; 
   else
      clusterNavFile= fullfile(pathname, filename);
      setappdata(hfig, 'clusterNavFile', clusterNavFile);
      btnSave_Callback(hObject, eventdata, handles);
   end   


% --- Executes on button press in checkboxLog.
function checkboxLog_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxLog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxLog
logScale = get(hObject,'Value');
hfig = handles.figMain;
hax = getappdata(hfig,'ldaAxes');
if logScale
  set(hax,'YScale','log')
else
  set(hax,'YScale','linear');
end


% --- Executes on button press in checkboxPairwiseUpdate.
function checkboxPairwiseUpdate_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxPairwiseUpdate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear_lda_plots(handles);
update_lda(handles);
% Hint: get(hObject,'Value') returns toggle state of checkboxPairwiseUpdate


