%%%% ASSIGN_FIDUCIALS: assign fiducial points between a reference and moving image with cpselect,
%%%% create a warped version of the moving image with projective warping, save the results
%
% assign_fiducials(movingimage,refimage,old_pair_file,savename,reverse_points)
%
%%%% last edited 20/5/5 on thermaltake

function assign_fiducials(movingimage,refimage,old_pair_file,savename,reverse_points)

transform_type = 'projective';
% tifpars.conv_factor = 16^4; % for savetif
tifpars.conv_factor = 1; % for savetif
tifpars.keepbitdepth = 1;
reload_refimage = 1; % if loading previous work, reload newest version of movingimage rather than stored version
reload_movingimage = 1; 
maintain_black_white = 1; % if movingimage contains only 2 pix vals (0 and something else), set output to have only these pix vals

[movingimage movingimageFile] = procImageInput('movingimage');
[refimage refimFile] = procImageInput('refimage');

if ~exist('savename','var') || isempty(savename)
    savename = [getfname(movingimageFile) '_reg_to_' getfname(refimFile)];
end

if exist([savename '.mat'],'file') %%% if this image pair has had points assigned
% % % % % % %     go_on = input(['Load and update file ' savename '? (y/n) '],'s');
    go_on = 'y';
    if strcmp(go_on,'y')
        r = refimage;
        m = movingimage;
        load([savename '.mat']); % load previous work on the same file pair
        if reload_refimage
            refimage = r;
        end
        if reload_movingimage
            movingimage = m;
        end
        clear r m
    else
        error('quitting function')
    end
    [movingimage refimage] = remove2ndPlane(movingimage,refimage);
    [movingpoints fixedpoints] = cpselect(movingimage,refimage,movingpoints,fixedpoints,'Wait',true); % start with preexisting points
elseif exist('old_pair_file','var') && ~isempty(old_pair_file) %%% if working from points found in a different file pair
    fprintf(['Loading points from ' getfname(old_pair_file) '...\n'])
    oldpoints = load(old_pair_file,'movingimage','refimage','movingpoints','fixedpoints');
    [movingimage refimage] = remove2ndPlane(movingimage,refimage);
    if ~exist('reverse_points','var') || ~reverse_points
        [movingpoints fixedpoints] = cpselect(movingimage,refimage,oldpoints.movingpoints,oldpoints.fixedpoints,'Wait',true); % start with preexisting points    
    elseif reverse_points % exchange old movingpoints and fixedpoints 
        [movingpoints fixedpoints] = cpselect(movingimage,refimage,oldpoints.fixedpoints,oldpoints.movingpoints,'Wait',true); % start with preexisting points    
    end
else
    [movingpoints fixedpoints] = cpselect(movingimage,refimage,'Wait',true);
end

geotransform = fitgeotrans(movingpoints,fixedpoints,transform_type);
worldref = imref2d(size(refimage));
movingimage_warped = imwarp(movingimage,geotransform,'OutputView',worldref);

if maintain_black_white
   unq_moving = unique(movingimage(:));
   if length(unq_moving) < 3 && unq_moving(1) == 0 % if original was black-white
       thresh = unq_moving(2);
       movingimage_warped(movingimage_warped >= thresh) = max(unq_moving); % binarize
       movingimage_warped(movingimage_warped < thresh) = unq_moving(1); % binarize
   end
end

timelastsaved = datestr(now);

save(savename)
savetif(movingimage_warped,savename,tifpars)

end



%%% remove extraneous dimensions to work with cpselect
function [im1_out, im2_out] =  remove2ndPlane(im1_in,im2_in) 
    if size(im1_in,3) ~= 3
        im1_out = im1_in(:,:,1);
    else
        im1_out = im1_in;
    end
    if size(im2_in,3) ~= 3
        im2_out = im2_in(:,:,1);
    else
        im2_out = im2_in;
    end
end