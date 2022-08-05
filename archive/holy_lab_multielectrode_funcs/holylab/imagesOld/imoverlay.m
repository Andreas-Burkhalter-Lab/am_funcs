function imout = imoverlay(im1,im2,color1,color2,options)
  if (nargin < 5)
    options.style = 'max';
  end
  [m,n] = size(im1);
  for i = 1:3
    imtmp1 = double(im1) * color1(i);
    imtmp2 = double(im2) * color2(i);
    if strcmp(options.style,'max')
      imout(:,:,i) = max(imtmp1,imtmp2);
    elseif strcmp(options.style,'sum')
      imout(:,:,i) = imtmp1+imtmp2;
    end
  end
  immax = max(imout(:));
  imout = imout / immax;