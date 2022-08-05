function output = lkInverseAffine(vid_file,startingPoints,startf,endf,needsDemosiac,s)
%% Lucas Kanade semi-live inverse compositional -- affine by John Boyle
% Code is used to find local deformations on a material.
%
% inputs: 
% vid_file - video file to load
% startingPoints - a n x 2 vector of center points [y , x] of positions of points to track
% startf - the first frame in the video we want to start tracking
% endf - the last frame in the video we want to track
% needsDemosiac - set true if video is bayer coded
% s - a name of the sample
% 
% outputs : (contained in a structure called 'output')
% points        - a cell array containing the location of the center of  
%                 each region in each frame 
% secondChances - an array which tells how many times each template was 
%                 updated with each value corresponding to each point in
%                 points
% F             - a cell array containing deformation gradient tensors for
%                 each point organized by F{frame}(1:9,point) = 
%                 [F11 F12 T1 F21 F22 T2 0 0 1]
%                                   
% Semi-live inverse compositional method: 
% By default, the code tries to match regions in the current video frame
% with the very first frame. By doing this, the code directly calculates the
% deformation of every region relative to the first frame. To aid in
% convergence, for each frame the deformation parameters from the previous
% frame are used as an initial guess. 
% Proceeding this way, however may lead to less convergence as the
% deformations become very large. To account for this the code allows the
% template to be updated dynamically in the case that a region fails to
% converge on a solution. It does this in the following steps:
%
% 1. If after max_iterations we still have: ||dP|| < epsilon , we consider 
%    the point as "failed" to converge.
% 2. When a point fails to converge, we then go back to the frame
%    immediately before and define a new template region. 
% 3. For the new template region we recalculate steepest descent images and
%    hessian matrices.
% 4. We also define a set of compositional parameters, which is the simply 
%    parameter values from the previous frame.    
% 5. For each failed point we then use the new template data to recalculate
%    the current warp parameters, which are now relative to the previous 
%    frame. 
% 6. If the region still fails to converge, we discard it (unTrusted)
% 7. When calculating the total deformation (F) relative to the first frame 
%    we compose the current parameters with the compositional
%    parameters then calculate the total deformation. 
% 8. For robustness, we track everytime a template has been updated (given
%    a secondChance). Each time the template is adjusted the region likely
%    contains more error. 
% 
% For simplicity, the current parameters are composed every frame with a 
% set of compositional parameters which are initially set to be 0. The
% compositional parameters are only updated if the optimization loop fails
% to converge. 
% 
% equation references from:
% Baker, S., & Matthews, I. (2004). Lucas-Kanade 20 Years On: A Unifying Framework. International Journal of Computer Vision, 56(3), 221-255. doi:10.1023/B:VISI.0000011205.11775.fd

% options
boxSize = [70 70]; % how big do we want to make the region tracking?
b = 5; % gradient padding (improves quality)

options.color = 2; % 4 - take mean, 1,2, or 3, use R G or B
options.fstep = 1; % frame step (change to only analyze every (fstep) frame)
options.use_par_cpu = false; % process points in parallel
options.max_iter = 80; % maximum number of optimization iterations
options.grad_sig = .2; % define the gradient step size (smaller takes longer to converge larger may overshoot)
options.plotLive = true; % set true to update a figure every loop with the current location of the tracked regions
options.verbose = true; % if true, the program will show a figure every optimization loop as it descends on the correct answer
options.error_threshold = 0.01; % stop iterating when: ||dP|| < error_threshold
options.nonTranslationIterations = 0; % enforce certain parameters over others in x number of the first iterations
options.secondChance = true; % if true, when a solution fails to converge the program will attempt to update the basis images, steepest descent images, and hessian matrices of the previously known location
options.waitEachLoop = false; % enable to trigger a pause every frame loop 
options.LevenMarq = false;

nframe = endf-startf; % number of frames to process

%% preprocessing the video 
movie = mmreader(vid_file); % create video object
image1 = read(movie,startf);
if needsDemosiac
    image1 = demosaic(image1(:,:,1),'gbrg'); % demosiac the image before displaying
end
if options.color == 4; % if color is 4 take mean of the image
    image1 = mean(image1,3); % mean 
    image1 = double(image1); % convert to double for processing
else % if not 4 then use the specific color
    image1= image1(:,:,options.color); % use specific color
    image1 = double(image1); % make double
end

frame = 1; % start index at first frame
points{1} = startingPoints; % create our points cell array
fprintf('sample: %s \t regions:\t %i \n',s,numel(points{1}(:,1)))

%% decide how to distribute the work between number of workers
% for parallel computation. if not using parallel, just use all of the
% points
if options.use_par_cpu
    clear ans % matlab has no nice way of getting number of current workers
    matlabpool size; % so do matlabpool size 
    poolsize=ans; %#ok<NOANS> % and get number of workers from the ans
    if poolsize==0; % if we have no workers
        matlabpool %start a new matlabpool
    end
    matlabpool size; % now get the size of the pool
    poolsize=ans;    %#ok<NOANS> % from the ans again
    numberperlab=length(points{frame})/poolsize; % divide our points between the workers
    % single program multiple data -  here we decide how to distribute the work between the
    % available workers
    spmd % now lets assign the points to each worker
        if labindex<poolsize % in this case input_points becomes a composite on each lab containing the points the lab should process
            input_points=1+(labindex-1)*fix(numberperlab):(labindex)*fix(numberperlab);
        else
            input_points=1+(labindex-1)*fix(numberperlab):length(points{frame});
        end
    end
else % if not using parallel cpu then just assign all points to a normal variable
    input_points = 1:length(points{frame});
end

%% precomputation
% precompute jacobian for Warp(x;0) 
% the jacobian will be the same for all images:
m = boxSize(1)+1 ; n = boxSize(2)+1; % m x n images from the boxsize
[xg,yg] = ndgrid(0:m-1,0:n-1); % create a grid of the points 
Tcenter = [ m n ]/2; % find the centerpoint
x = xg-Tcenter(1); % set the centerpoint equal to 0 
y = yg-Tcenter(2); % " " 

xn = numel(x); % how many pixels do we have 
O = zeros(size(x(:))); % zero fill for jacobian
l =  ones(size(x(:))); % one fill for jacobian 
% building jacobian for affine warp at W(x;0) (equation 8) 
jac_x = [x(:) O    y(:) O    l O ]; % build jacobian for x
jac_y = [O    x(:) O    y(:) O l ]; % build jacobian for y

% precomputation of the hessian matrices and the steepest descent images
% these are unique for each region so lets loop through them
% also will be prebuilding matrices that will be used later. 
if options.use_par_cpu
    spmd % on each worker:
        points1 = points{1}(input_points,:); % get the starting location of the points for tracking
        % prebuild some matrices on each worker for use later
        Pc = zeros(6,length(input_points)); % prebuild Parameter matrix
        Fc = zeros(9,length(input_points)); % prebuild deformation gradient tensor matrix
        i_count_c = zeros(length(input_points),1); % prebuild a count of iterations
        % define each template box region for tracking
        cpbox(:,:) = [points1(:,:)-repmat(floor(boxSize/2),length(input_points),1) repmat(boxSize,length(input_points),1)]; %#ok<UNRCH>
        % prebuild greatest descent, hessian, and diagonal hessian (Hlm) matrices
        H = zeros(6,6,numel(input_points)); % hessian
        Hlm = zeros(6,6,numel(input_points)); % diagonal hessian
        dTdWdP = zeros(xn,6,numel(input_points)); % steepest descent 
        dTdWdP_ = zeros(xn,6); % temporary steepest descent
        bIms = zeros(m,n,numel(input_points)); % create a 3d matrix where we will store all of the template images
        Pcomp = zeros(6,numel(input_points)); % build an array of compositional P 
        % loop through each region and precompute the steepest descent,
        % hessian, and diagonal hessian
        for i = 1:numel(input_points)
        % precompute steepest descent and inverse hessian
        [dTdWdP(:,:,i),H(:,:,i),Hlm(:,:,i)] = LKpreComputeLM( ...
            image1(cpbox(i,2)-b:cpbox(i,2)+cpbox(i,4)+b,cpbox(i,1)-b:cpbox(i,1)+cpbox(i,3)+b) ...
            ,b,xn,jac_x,jac_y,dTdWdP_,options); 
        % save base images in 3d array
        bIms(:,:,i) = image1(cpbox(i,2):cpbox(i,2)+cpbox(i,4),cpbox(i,1):cpbox(i,1)+cpbox(i,3));
        end
    end
else % if not using parrallel cpu (code should be identical to that in spmd)
        points1 = points{1}(input_points,:); % get the starting location of the points for tracking
        % prebuild some matrices on each worker for use later
        Pc = zeros(6,length(input_points)); % prebuild Parameter matrix
        Fc = zeros(9,length(input_points)); % prebuild deformation gradient tensor matrix
        i_count_c = zeros(length(input_points),1); % prebuild a count of iterations
        % define each template box region for tracking
        cpbox(:,:) = [points1(:,:)-repmat(floor(boxSize/2),length(input_points),1) repmat(boxSize,length(input_points),1)]; %#ok<UNRCH>
        % prebuild greatest descent, hessian, and diagonal hessian (Hlm) matrices
        H = zeros(6,6,numel(input_points)); % hessian
        Hlm = zeros(6,6,numel(input_points)); % diagonal hessian
        dTdWdP = zeros(xn,6,numel(input_points)); % steepest descent 
        dTdWdP_ = zeros(xn,6); % temporary steepest descent
        bIms = zeros(m,n,numel(input_points)); % create a 3d matrix where we will store all of the template images
        Pcomp = zeros(6,numel(input_points)); % build an array of compositional P 
        % loop through each region and precompute the steepest descent,
        % hessian, and diagonal hessian
        for i = 1:numel(input_points)
        % precompute steepest descent and inverse hessian
        [dTdWdP(:,:,i),H(:,:,i),Hlm(:,:,i)] = LKpreComputeLM( ...
            image1(cpbox(i,2)-b:cpbox(i,2)+cpbox(i,4)+b,cpbox(i,1)-b:cpbox(i,1)+cpbox(i,3)+b) ...
            ,b,xn,jac_x,jac_y,dTdWdP_,options); 
        % save base images in 3d array
        bIms(:,:,i) = image1(cpbox(i,2):cpbox(i,2)+cpbox(i,4),cpbox(i,1):cpbox(i,1)+cpbox(i,3));
        end
end

% if we are showing an updated figure every loop, we need to show the 
% positions of the original points
if options.plotLive
   figure(1)
   close all
   current_im =imagesc(image1);
   set(gca,'NextPlot','replacechildren') 
   hold on
   current_plot = plot(gca,points{1}(:,1),points{1}(:,2),'b.');
   colormap('gray')
   set(gca,'NextPlot','replacechildren')
end

%% frame processing loop
bf = startf; % current search image
bframe = 1; % base frame
cframe = 2; % current frame
trusted = 1:length(points{frame}); % prebuild an array of what we trust (everything)
secondChances = zeros(1,size(points{frame})); % prebuild an array of which points have had their templates recalculated
untrusted = []; % nothing is untrusted yet, create an empty array
b_frac_offset = zeros(size(points{frame}));
% cframe == current frame (relative to 1st frame of analysis)
% bframe == base frame (relative to 1st frame of analysis)
% cf == current frame in the video (absolute frame)
% bf == base frame in the video (absolute)
% base_im == loaded base image
% curr_im == loaded current image
% set P initially to all zeros
Pa{1} = zeros(6,size(points{1},1));
while cframe < nframe+10 % number of frames
    tic % how long does it take to do the frame?
    cf = startf + (cframe-1)*options.fstep ; % current frame (relative to
    % setup current frame)
    if cframe ~=2 % check to see if we can just use cf for bf
        last_im = curr_im; % instead of reloading the base image, just use the one from before
    end    
    curr_im = read(movie,cf); % current frame
    % demosaic, converto to greyscale, and make the images double
    if needsDemosiac % is our video bayer coded?
        curr_im = demosaic(curr_im(:,:,1),'gbrg');
    end
    if options.color == 4; % 4 means use the mean of the image
        curr_im = mean(curr_im,3);
        curr_im = uint8(curr_im);
        curr_im = double(curr_im);
    else % if not 4, use the explicitly defined color
        curr_im = curr_im(:,:,options.color);
        curr_im = double(curr_im);
    end
    fprintf('frame: %i \t bframe: %i \t cf: %i \t bf: %i \n',cframe, bframe, cf, bf)
    % start the loop of tracking each point
    if options.use_par_cpu
        spmd
            % make the current boxes
            cpboxcurr(:,:) = [round(points{cframe-1}(input_points,:))-repmat(floor(boxSize/2),length(input_points),1) repmat(boxSize,length(input_points),1)]; %#ok<UNRCH>
            % calculate the change between images.
            [Pc Fc i_count_c] = LKloopLM(bIms,curr_im,dTdWdP,H,Hlm,1:numel(input_points),Pc,Fc,Pcomp,i_count_c,cpboxcurr,xg,yg,options,cframe);            
            i_count_nsc = i_count_c; % save the values before updating to identify which points have been updated, nsc: no second chance
            if cframe>2 && options.secondChance % only correct if we are not in the first frame
                % find which failed to track
                failed = find(i_count_c == options.max_iter)';
                % lets try and update them using the last steepest descent image
                % to see if we can recover them
                for i = failed
                    % precalculation loop, using last im in the location from last
                    % frame
                    [dTdWdP(:,:,i),H(:,:,i),Hlm(:,:,i)] = LKpreComputeLM( ...
                        last_im(cpboxcurr(i,2)-b:cpboxcurr(i,2)+cpboxcurr(i,4)+b,cpboxcurr(i,1)-b:cpboxcurr(i,1)+cpboxcurr(i,3)+b) ...
                        ,b,xn,jac_x,jac_y,dTdWdP_,options); % precompute steepest descent and inverse hessian
                    % save the base images in a 3d array
                    bIms(:,:,i) = last_im(cpboxcurr(i,2):cpboxcurr(i,2)+cpboxcurr(i,4),cpboxcurr(i,1):cpboxcurr(i,1)+cpboxcurr(i,3));
                    Pc(:,i) = zeros(6,1); % update Pc so that the new warp starts fresh!
                    Pcomp(:,i) = Pprev(:,i); % update our composite with the previous confident warp
                    Pcomp([5 6],i) = 0; % translation of the composite is always 0!
                end
                % now update the variables
                [Pc Fc i_count_c] = LKloopLM(bIms,curr_im,dTdWdP,H,Hlm,failed,Pc,Fc,Pcomp,i_count_c,cpboxcurr,xg,yg,options,cframe);
                %failed = find(i_count_c == options.max_iter)';
            end
            Pprev = Pc; % save Pc for potential updates next loop 
        end
        % collect data from workers
        for j=1:numel(Pc)
            lengthx(j,:)=size(Pc{j});
        end
        bis=0;
        for j=1:numel(Pc)
            von=bis+1;
            bis=bis+lengthx(j,2);
            i_count_nsc_temp(:,von:bis) = i_count_nsc{j}'; 
            i_count_temp(:,von:bis) = i_count_c{j}';
            Ptemp(:,von:bis) = Pc{j};
            Ftemp(:,von:bis) = Fc{j};
        end
        
        % assign variables
        icountnsc{cframe} = i_count_nsc_temp;
        Pa{cframe} = Ptemp;
        Fa{cframe} = Ftemp;
        icount{cframe} = i_count_temp;
    else % not using parallel CPU
        % make the current boxes
        cpboxcurr(:,:) = [round(points{cframe-1}(input_points,:))-repmat(floor(boxSize/2),length(input_points),1) repmat(boxSize,length(input_points),1)]; %#ok<UNRCH>
        % calculate the change between images.
        [Pc Fc i_count_c] = LKloopLM(bIms,curr_im,dTdWdP,H,Hlm,1:numel(input_points),Pc,Fc,Pcomp,i_count_c,cpboxcurr,xg,yg,options,cframe);
        i_count_nsc = i_count_c; % save the values before updating to identify which points have been updated, nsc: no second chance      
        if cframe>2  && options.secondChance % only correct if we are not in the first frame
            % find which failed to track
            failed = find(i_count_c == options.max_iter)';
            % lets try and update them using the last steepest descent image
            % to see if we can recover them
            for i = failed
                % precalculation loop, using last im in the location from last
                % frame
                [dTdWdP(:,:,i),H(:,:,i),Hlm(:,:,i)] = LKpreComputeLM( ...
                    last_im(cpboxcurr(i,2)-b:cpboxcurr(i,2)+cpboxcurr(i,4)+b,cpboxcurr(i,1)-b:cpboxcurr(i,1)+cpboxcurr(i,3)+b) ...
                    ,b,xn,jac_x,jac_y,dTdWdP_,options); % precompute steepest descent and inverse hessian                
                % save the base images in a 3d array
                bIms(:,:,i) = last_im(cpboxcurr(i,2):cpboxcurr(i,2)+cpboxcurr(i,4),cpboxcurr(i,1):cpboxcurr(i,1)+cpboxcurr(i,3));
                Pc(:,i) = zeros(6,1); % update Pc so that the new warp starts fresh!
                Pcomp(:,i) = Pprev(:,i); % update our composite with the previous confident warp
                Pcomp([5 6],i) = 0;
            end
            % now update the variables
            [Pc Fc i_count_c] = LKloopLM(bIms,curr_im,dTdWdP,H,Hlm,failed,Pc,Fc,Pcomp,i_count_c,cpboxcurr,xg,yg,options,cframe);

        end
        Pprev = Pc; % save Pc for potential updates next loop 
        
        % for consistency... like collecting data from workers
        Ptemp(:,:) =  Pc;
        Ftemp(:,:) = Fc;
        i_count_temp(:,:) = i_count_c;
        i_count_nsc_temp(:,:) = i_count_nsc'; 
        
        
        % assign variables
        icount{cframe} = i_count_temp;
        Pa{cframe} = Ptemp;
        Fa{cframe} = Ftemp;
        icountnsc{cframe} = i_count_nsc_temp;
    end
    
    % find which we couldn't track
    frameSecondChance = icountnsc{cframe}==options.max_iter;
    secondChances = secondChances + frameSecondChance;
    framefailed = find(icount{cframe}==options.max_iter);        
    untrusted = unique([framefailed(:) ; untrusted(:)]); % add anything that failed to converge
    trusted = setdiff(trusted,untrusted); % pop them frm what we believe 
    recovered = intersect(trusted,find(frameSecondChance==1));
    t_toc(cframe) = toc;
    fprintf('fr: %d \t elapsed: %3.2f \n',cframe,t_toc(cframe))
    in_frac_offset = points{cframe-1}-round(points{cframe-1});
    b_frac_offset(framefailed,:) = points{cframe-1}(framefailed,:)-round(points{cframe-1}(framefailed,:));
    points{cframe} = points{cframe-1}-Fa{cframe}([8 7],:)' - in_frac_offset - b_frac_offset;
    
    
    if options.plotLive
        figure(1)
        current_im = imagesc(curr_im);
        hold on
        new_plot = plot(gca,points{cframe-1}(trusted,1),points{cframe-1}(trusted,2),'b.');
        new_plot2 = plot(gca,points{cframe-1}(untrusted,1),points{cframe-1}(untrusted,2),'r.');
        now_plot3 = plot(gca,points{cframe-1}(recovered,1),points{cframe-1}(recovered,2),'m.');
        clear current_plot
        colormap gray
        drawnow;
        if options.waitEachLoop; pause; end
    end
    cframe = cframe + 1;   
end

pTrust = 1-length(untrusted)/(length(untrusted)+length(trusted)); % output the location of each point

%% saving the values for later
output.points = points; % save the locations of all the points to our output
output.secondChances = secondChances; % save how many times for each point the templates were recalculated
output.F = Fa; % save the deformation gradient tensors
