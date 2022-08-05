%%%% proc_image_inputs
% input is the name of an image arg to fill from caller workspace; can be
% filename or raw image; if empty, use uigetfile to select file
% [imageOut imageFilename] = procImageInput(queryName)
% updated 3/3/18 on thermaltake

function [imageOut imageFilename] = procImageInput(queryName)

if evalin('caller',['exist(''',queryName,''',''var'')'])... % if image arg exists in caller workspace
    && evalin('caller',['~isempty(',queryName,')']) %%% and is not empty
    imageIn = evalin('caller',queryName); % get image arg from caller workspace
else % if image isn't specified in caller yet
    imageIn = uigetfile('*',['Select file for input ''' queryName '''']); % select an image file
end

if ischar(imageIn) % if imageIn is a filename
    if exist(imageIn,'file')
        imageFilename = imageIn;
        imageOut = imread(imageIn);
    else 
        imageFilename = uigetfile('*',['File ' imageIn ' not found. Select file for argument ''' queryName '''.']);
        imageOut = imread(imageIn);
    end
elseif isnumeric(imageIn) % if imageIn is the actual image
    imageFilename = '';
    imageOut = imageIn;
else
    error(['Input arg ''' queryName ''' is not an image or filename.'])
end
        