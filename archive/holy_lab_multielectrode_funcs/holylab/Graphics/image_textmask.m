function textmask = image_textmask(sz,str,varargin)
  % image_textmask: create a mask of defined size using text
  %
  % This function is useful if you want to add text to pixels of an image
  %
  % Syntax:
  %   textmask = image_textmask(sz,str);
  %   textmask = image_textmask(sz,str,'FontSize',12,'position',[30 1],'VerticalAlignment','top',...);
  % Here,
  %   sz is a vector giving the number of pixels along the
  %     [vertical horizontal] axes (the image size), or in the general case
  %     the size of the final array.
  %   str is the text you want to display
  %   You can include property/value pairs to control placement, etc.
  % On output
  %   textmask is a logical array of size sz, true in the elements that
  %     correspond to the pixels of the text.
  %
  % You can see the text using
  %   imagesc(textmask)
  % or apply it with
  %   myimage(textmask) = uint8(255)
  % if myimage is a uint8 image.

  % This was slightly modified from the following posting:
  %    http://www.mathworks.com/support/solutions/en/data/1-1BALJ/index.html?solution=1-1BALJ
  
  pv = {'position',[1 sz(1)],'VerticalAlignment','top'};  % default parameter-value pairs
  if nargin > 2
    pv = [pv varargin];
  end
  
  % Create the figure with an axis of the correct size
  hf = figure('color','white','units','normalized','position',[0 0 1 1]);
  image(ones(sz))
  set(gca,'units','pixels','position',[5 5 sz(2)-1 sz(1)-1],'visible','off')

  % Draw the text
  text('units','pixels','string',str,pv{:});

  % Capture the text image
  % Note that the size will have changed by about 1 pixel
  pause(1);  % give enough time to render the figure
  textim = getframe(gca);
  textmask = textim.cdata(:,:,3) == 0;
  if length(sz) > 2
    textmask = repmat(textmask,[1 1 sz(3:end)]);
  end
  
  % Close the image
  close(hf)
end
