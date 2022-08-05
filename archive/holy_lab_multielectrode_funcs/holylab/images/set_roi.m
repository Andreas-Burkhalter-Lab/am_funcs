function roi_defs=set_roi(fig, time, roi_defs)
   rois_alltime=getappdata(fig, 'rois');
   if(isempty(rois_alltime))
      rois_alltime.time     =time;
      rois_alltime.tform    =roi_defs.tform;
      rois_alltime.defs_orig=roi_defs;
   else 
      if(isfield(roi_defs, 'tform_changed') && roi_defs.tform_changed)
         idxToReplace=find_index_to_replace(rois_alltime.time, time);
         if(idxToReplace<0) % if, insert in fact
            idxToIns=-idxToReplace;  % Here's where it will be inserted
            % We're going to copy, and we can't copy beyond the end
            idxToBreak = min([idxToIns length(rois_alltime.time)]);
            rois_alltime.time=rois_alltime.time([1:idxToBreak idxToBreak:end]);
            rois_alltime.tform=rois_alltime.tform([1:idxToBreak idxToBreak:end]);
            idxToReplace=idxToIns;
         end
         rois_alltime.time(idxToReplace)=time;
         rois_alltime.tform(idxToReplace)=roi_defs.tform;
      end
      
      % here assume the .action is consistent with its data, i.e. 'none' is
      % no change so it is fine to overwrite.
      labelRoiToDel=[];
      for idxRoi=1:length(roi_defs.label)
         if(strcmp(roi_defs.action{idxRoi}, 'deleted'))
            labelRoiToDel=[labelRoiToDel roi_defs.label(idxRoi)];
         end
      end
      roi_defs=delete_roi_defs(roi_defs, labelRoiToDel);
      rois_alltime.defs_orig=roi_defs;
   end
   
   % remove 4 fields in defs_orig (.tform_changed, .tform, .action, .handle )
   rois_alltime.defs_orig=rmfield_g(rois_alltime.defs_orig, {'tform', 'tform_changed', 'action','handle'} );
   
   % change .tform_changed and .action in returned roi_defs 
   roi_defs.tform_changed=0;
   for idx=1:length(roi_defs.label)
      roi_defs.action{idx}='none';
   end
   
   setappdata(fig, 'rois', rois_alltime);
   
   
% todo: separate this func   
function s=rmfield_g(s, fieldsToDel)
   idxValidField=[];
   for idx=1:length(fieldsToDel)
      if(isfield(s, fieldsToDel{idx}))
         idxValidField=[idxValidField idx];
      end
   end
   s=rmfield(s, fieldsToDel(idxValidField));

% todo: separate this func   
% return negative value if insert instead 
% alltime must be sorted.
function idxToReplace=find_index_to_replace(alltime, time)
   if(isempty(alltime))
      idxToReplace=1; 
      return;
   end
   
   idxToReplace=find(alltime==time);
   if(length(idxToReplace)>1)
      error('code error in set_roi.m');
   elseif(length(idxToReplace)==1)
      return;
   end
   
   if(time<alltime(1))
      idxToReplace=-1;
      return;
   end
   
   if(time>alltime(end))
      idxToReplace=-(length(alltime)+1);
      return;
   end
   
   idxSmaller=find(alltime<time);
   idxToReplace=-(idxSmaller(end)+1);

