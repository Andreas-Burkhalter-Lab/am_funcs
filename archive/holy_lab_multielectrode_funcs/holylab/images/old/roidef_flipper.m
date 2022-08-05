function roidef_flipper(roidefilename)
%roidef_flipper takes roi definitions files and flips their y axis
%orientation. This can be done just for fun or, more usefully, if you used
%stack_roi during the month of September, 2006 when we had a bug that
%transposed the vertical roi matrix orientation to that of the image data.
%The function automatically saves a new version of the roidef file with an
%altered name: [baseroideffilename '_flip.roidef].
%
%syntax: roidef_flipper(roidefilename)
%
%Copywrite 2006 Terrence Holekamp

load(roidefilename, '-mat')

height = header.height;

for i = 1:length(roi_defs)
    oldposition = roi_defs(i).posInPixel(2);
    newposition = (height - oldposition) + 1;
    roi_defs(i).posInPixel(2) = newposition;
end

[pathstr, filename] = fileparts(roidefilename);

save([filename '_flip.roidef'],'header','roi_defs', 'pixelPerUm', 'tform_info', '-mat');
    