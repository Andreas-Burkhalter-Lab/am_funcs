function roi_defs1=merge_roi_defs(roi_defs1, roi_defs2)
   if(isempty(roi_defs1))
      roi_defs1=roi_defs2;
      return;
   end
   if(isempty(roi_defs2))
      return;
   end
   
   names=fieldnames(roi_defs1);
   for idx=1:length(names)
      name=names{idx};
      switch name
         case get_roi_shared_fields % such as: tform
            % do nothing
         otherwise
            roi_defs1.(name)=[roi_defs1.(name) roi_defs2.(name)];
      end
   end

