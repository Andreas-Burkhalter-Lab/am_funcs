%%%% display movie while display simultaneous data trace
% in data trace, time must proceed down rows
%%% updated 2019-6-4 on thermaltake

function show_movie_with_timecourse(image_stack, data_trace, frames_per_sec)

save_into_movie = 1;
video_savename = 'saved_video';
data_plot_position = [0.1 0.05 0.8 0.35];
movie_plot_position = [0.2 0.5 0.6 0.35]; 

if iscell(image_stack) % if cell array, should be 1-D cell array where each cell contains one frame; frames should be equal size
    image_stack = reshape(cell2mat(reshape(image_stack,1,length(image_stack))), [size(image_stack{1}), length(image_stack)]); % turn into 3D matrix
end
    
if size(image_stack,3) ~= length(data_trace)
    error('movie and data trace must be same length')
end
nframes = size(image_stack,3);

frames_per_sec = vardefault('frames_per_sec',0.3);
ff = figure('units','normalized','outerposition',[0.2 .1 .6 .9 ]);
movie_plot = subplot('Position',movie_plot_position);
data_plot = subplot('Position', data_plot_position);
mindata = min(data_trace(:));
maxdata = max(data_trace(:));

frameperiod = 1/frames_per_sec;
for iframe = 1:nframes
    subplot(movie_plot)
    imagesc(image_stack(:,:,iframe))
    colormap(gca,gray)
    subplot(data_plot)
    plot(data_trace(1:iframe,:))
    xlim([0 nframes])
    ylim([mindata maxdata])
    ylabel('Locomotion (m/s)')
    title(['frame ' num2str(iframe) '/' num2str(nframes)])
    if save_into_movie
        framesave(iframe) = getframe(gcf);
    else
        pause(frameperiod)
    end
end

if save_into_movie
    writerObj = VideoWriter(video_savename);
    writerObj.FrameRate = frames_per_sec;
    open(writerObj);
    for i=1:length(framesave)
        thisframe = framesave(i);
        writeVideo(writerObj, thisframe);
    end
    close(writerObj);
end