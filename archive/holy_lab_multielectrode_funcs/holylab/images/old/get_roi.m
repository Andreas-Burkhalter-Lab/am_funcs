function rois=get_roi(fig, time)
   rois_alltime=getappdata(fig, 'rois');
   if(isempty(rois_alltime))
      rois=[];
      return;
   end
   
   rois=rois_alltime.defs_orig;
   
   % now need add .tform field
   
   alltimes=rois_alltime.time;
   
   % if not in range:
   if(time<=alltimes(1))
      rois.tform=rois_alltime.tform(1);
      return
   elseif(time>=alltimes(end))
      rois.tform=rois_alltime.tform(end);
      return;
   end
   
   % if exactly match:
   idx=find(alltimes==time);
   if(~isempty(idx))
      rois.tform=rois_alltime.tform(idx);
      return;
   end
   
   % if in range but not exactly match
   idx=find(alltimes<time);
   prevTime=alltimes(idx(end));
   nextTime=alltimes(idx(end)+1);
   prevTform=rois_alltime.tform(idx(end));
   nextTform=rois_alltime.tform(idx(end)+1);
   rois.tform=interpolate_tform([prevTime nextTime], [prevTform nextTform], time);
   
   
function T=interpolate_tform(timeFromAndTo, tforms, time)
   translation1=tforms(1).tdata.T(end, 1:2);
   translation2=tforms(2).tdata.T(end, 1:2);
   translation=translation1+(translation2-translation1)*(time-timeFromAndTo(1))/diff(timeFromAndTo);
   affine_matrix=[eye(2); translation];
   T=maketform('affine', affine_matrix);
   
   





