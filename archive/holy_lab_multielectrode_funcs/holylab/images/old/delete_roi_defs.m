function roi_defs=delete_roi_defs(roi_defs, labelsToDel)
   [labelToKeep, idxToKeep]=setdiff(roi_defs.label, labelsToDel);
   
   names=fieldnames(roi_defs);
   for idx=1:length(names)
      name=names{idx};
      switch name
         case get_roi_shared_fields % 
            % do nothing
         otherwise
            roi_defs.(name)=roi_defs.(name)(idxToKeep);
      end
   end
