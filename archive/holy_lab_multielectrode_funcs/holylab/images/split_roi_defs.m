function cellarray_roi_defs=split_roi_defs(roi_defs)
% SPLIT_ROI_DEFS: convert ROI structure to cell array of single ROIs
% Syntax:
%   cellarray_roi_defs=split_roi_defs(roi_defs)
% where
%   roi_defs has the format defined in ROISTRUCT;
%   cellarray_roi_defs is a cell array of length nrois, where each
%     element has the format defined in ROISTRUCT, but with only a single
%     roi.
  
   cellarray_roi_defs={};
   if(isempty(roi_defs))
      return;
   end
   
   names=fieldnames(roi_defs);
   for idxRoi=1:length(roi_defs.label)
      for idx=1:length(names)
         name=names{idx};
         switch name
            case get_roi_shared_fields % such as: tform
               cellarray_roi_defs{idxRoi}.(name)=roi_defs.(name);
            otherwise
               cellarray_roi_defs{idxRoi}.(name)=roi_defs.(name)(idxRoi);
         end
      end
   end
   
