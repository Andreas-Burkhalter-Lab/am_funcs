%% how do we get a target frame? 

% T = meanImg; % this is cheating. Using the known average frame
% T = mean(mov,3); % how about using the average frame
% T = mov(:,:,20); % how about using a single frame? 
% T = findTarget(mov); % how about using an average of a few highly correlated frames
T = findTargetIterative(mov); % how about iteratively aligning a subset of frames

% determine y and x offsets
nFramesPerBatch = 1000;
nBatches = ceil(size(mov,3)/nFramesPerBatch);


% run the cross-correlation 
clear xoffset yoffset
for j = 1:nBatches
    trange = (j-1)*nFramesPerBatch + [1:nFramesPerBatch];
    trange(trange>size(mov,3)) = [];
    
%     [yoffset(trange), xoffset(trange)] = crossCorrelation(T, mov(:,:,trange));
    [yoffset(trange), xoffset(trange)] = ...
        crossCorrelation(T, mov(:,:,trange), 'phase');
end

% do these match the known ground truth offsets?
c1 = corrcoef(yoffset, xyshift(1:numel(yoffset),1));
c2 = corrcoef(xoffset, xyshift(1:numel(yoffset),2));

fprintf('Correlation between computed and true offsets\n')
fprintf('Y: %2.2f, X:%2.2f \n', c1(1,2), c2(1,2))

figure('Position', [100 100 800 400])
plot(xyshift(:,1), 'Linewidth', 2)
hold all
plot(yoffset)
hold all

legend('true', 'computed')

%% shift frames by the determined offsets
for j = 1:size(mov,3)
    % shift each frame *back* by the determined amounts
    mov(:,:,j) = shift_frame(mov(:,:,j), -yoffset(j), -xoffset(j));
end

%% check that registration worked
% pause('Press a key+enter to show registered movie')

figure('Position', [100 100 800 700]); 
colormap('gray')

for j = 1:200
    imagesc(mov(:,:,j), [0 4000])
    title('registered movie')
    drawnow
    pause(0.05)
end

%%