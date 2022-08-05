function rois=parse_rois(roifile)
% parse roi file saved by roidef

   [status, txt]=load_text_file(roifile);
   if(status~=0)
      error('failed to open roi def file');
   end
   
   lines=cellstr(split_str(txt, char(10)));
   for roiIndex=1:length(lines)
      line=lines{roiIndex};
      line=strrep(line, ';', char(10));
      rois(roiIndex).ver=key2value(line, 'roidef');
      rois(roiIndex).label=key2value(line, 'label');
      rois(roiIndex).from=str2num(key2value(line, 'from'));
      rois(roiIndex).to=str2num(key2value(line, 'to'));
      rois(roiIndex).active=str2num(key2value(line, 'active'));
      rois(roiIndex).geo=str2num(key2value(line, 'geo'));
   end
