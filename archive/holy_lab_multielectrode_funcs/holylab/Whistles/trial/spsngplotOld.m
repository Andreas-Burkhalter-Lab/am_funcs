function spsngplot(varargin)
% SPSNGPLOT: display a sparse sonogram (uses decimation)
% Usage:
%   spsngplot(sng)
%   spsngplot(sng,f,t)   (to specify frequency and time axes)
  if isnumeric(varargin{1})
    action = 'setup';   % The first call
  else
    action = varargin{1};
    hax = varargin{2};
  end
  switch action
   case 'viewrange',
    newrange = varargin{3};
    data = get(hax,'UserData');
    totalrange = size(data.sng,2);
    if any(newrange > totalrange)
      error('Desired range exceeds total size')
    end
    axpos = get(hax,'Position');
    npix = axpos(4);
    sngd = data.sng(:,newrange);
    if (diff(newrange) > npix)
      sngd = spsngdec(sngd,npix);
    end
    imagesc(data.f,data.t(newrange),log10(abs(sngd)));
    axis xy;
    colormap(1-gray);
   case 'setup',
    sng = varargin{1};
    if (nargin > 1)
      f = varargin{2};
    else
      f = 1:size(sng,1);
    end
    if (nargin > 2)
      t = varargin{2};
      if (length(t) < size(sng,2))
        t = linspace(min(t),max(t),size(sng,2));
      end
    else
      t = 1:size(sng,2);
    end
    data.sng = sng;
    data.f = f;
    data.t = t;
    set(gca,'UserData',data,'NextPlot','replacechildren');
    spsngplot('viewrange',gca,1:size(sng,2));
  end
