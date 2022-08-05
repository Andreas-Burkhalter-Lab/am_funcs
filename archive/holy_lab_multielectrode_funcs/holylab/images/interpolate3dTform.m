function tform=interpolate3dTform(tform_info, idxStack)
% NOTE: idxStack is 0-based
   if(isempty(tform_info.tform_def_times))
      tform=unit_tform(3);
      return;
   end % if, no tform was defined
   
   % add the unit tform
   tform_info.tform_def_times(end+1)=tform_info.roi_def_time;
   tform_info.tforms=[tform_info.tforms unit_tform(3)];
   
   % now sort by time
   [tform_info.tform_def_times, indices]=sort(tform_info.tform_def_times);
   tform_info.tforms=tform_info.tforms(indices);
   
   
   if(idxStack<=tform_info.tform_def_times(1))
      tform=tform_info.tforms(1);
      return;
   end
   
   if(idxStack>=tform_info.tform_def_times(end))
      tform=tform_info.tforms(end);
      return;
   end

   idx=find(tform_info.tform_def_times==idxStack);
   if(~isempty(idx))
      tform=tform_info.tforms(idx);
      return;
   end

   % TODO: next code can be simplified b/c time is sorted
   
   % now come to the real interpolation
   indices=find(tform_info.tform_def_times<idxStack);
   timeLeft=max(tform_info.tform_def_times(indices));
   indices=find(tform_info.tform_def_times>idxStack);
   timeRight=min(tform_info.tform_def_times(indices));
   
   idxLeft=find(tform_info.tform_def_times==timeLeft);
   idxRight=find(tform_info.tform_def_times==timeRight);
   
   tformLeft =tform_info.tforms(idxLeft);
   tformRight=tform_info.tforms(idxRight);
   
   curTime=idxStack;
   
   [A1, B1]=revert_maketform(tformLeft);
   [A2, B2]=revert_maketform(tformRight);
   
   factor=(curTime-timeLeft)/(timeRight-timeLeft);
   
   A=factor*(A2-A1)+A1;
   B=factor*(B2-B1)+B1;
   
   tform=maketform('affine',[A; B]);

   
function [A,B]=revert_maketform(tform)
   data=tform.tdata.T;
   A=data(1:3,1:3);
   B=data(4,  1:3);
   
