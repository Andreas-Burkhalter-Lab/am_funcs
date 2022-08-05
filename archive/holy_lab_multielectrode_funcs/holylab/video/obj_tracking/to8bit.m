function result=to8bit(im)
   im=double(im);
   minValue=min(im(:));
   result=(im-minValue)/(max(im(:))-minValue)*255;
   