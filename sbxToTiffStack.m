 %%% write all sbx files in current directory in order of filename into a single tiff stack
 %%%    or enter second arg to specify .sbx file list to write into stack
 % updated 2/12/18 
 
 function [] = sbxToTiffStack(outputname,dd)
 
 
 if ~exist('dd','var') || isempty(dd)
     dd = dir;
     dd= struct2table(dd);
    dd = dd.name(:,:);
    dd = dd(3:end); % remove top blank entries
    dd(cellfun(@(x)length(x)<4,dd)) = [];
    issbx = cellfun(@(x)strcmp(x(end-3:end),'.sbx'),dd,'UniformOutput',0);
    issbx = cell2mat(issbx);
    dd = dd(issbx);
 end
    
e = squeeze(sbxread(getfname(dd{1}),0,1));
imwrite(e,[outputname '.tif']);

for i = 2:length(dd)
    e = squeeze(sbxread(getfname(dd{i}),0,1));
    imwrite(e,[outputname '.tif'],'WriteMode','append');
end