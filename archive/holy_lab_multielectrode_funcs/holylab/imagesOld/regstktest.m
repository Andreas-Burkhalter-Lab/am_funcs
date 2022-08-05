function [pos,mov] = regstktest(img,rect,nframes,mag,pos)
  % Given a single frame img, a "window" rect defining a subset of the
  % image, create a new movie with this image jiggled around under the
  % window.
  % Syntax on rect is [left bottom width height]
  if (nargin < 5)
    dpos = mag*randn(nframes-1,2);
    pos = round([0 0; cumsum(dpos,1)]);
  end
  if (nargout < 2)
    return;   % Fast processing if user just wants to see amount of drift
  end
  for i = 1:nframes
    % Left column is x, right column is y
    rng = [rect(1:2); rect(1:2)+rect(3:4)] + [pos(i,:);pos(i,:)];
    mov(i).data = img(rng(1,2):rng(2,2),rng(1,1):rng(2,1));
  end
  