function image = flattentifstack()
% flattentifstack opens a gui-led process to flatten 16-bit tif stacks
%
% syntax: flattentifstack ; from any line
% See also 

% Copyright 2008 Julian P. Meeks

basefile = UIGetFiles('*.tif', 'Please select files to be flattened', pwd);
if size(basefile,2) < 1
    error('No input file selected, exiting program');
    image = [];
    return;
end

for idx = 1:size(basefile,2)
    temp = tiffread(basefile{idx});
    % (IN THEORY, checks for consistency go here)
    img(:,:,idx) = temp.data;
end

% TO DO: give UI options choice for type of flattening to do

image = max(img,[],3);

figure('position', [10 10 0.75*1024 0.75*1024]);
imagesc(image);
colormap(gray);

answer = questdlg('Does this look right?', 'Verify image');
if answer == 'Yes'
    imwrite(image,[basefile{1}(1:end-4) '.flat.tif'],  'tif');
end

end


