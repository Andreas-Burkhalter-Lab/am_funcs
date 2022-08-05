function discrep = regstkchk(mov,pos)
  nframes = length(mov);
  imref = mov(1).data;
  sz = size(imref);
  xindx = 1:sz(2);
  yindx = 1:sz(1);
  for i = 1:nframes
    xindxn = xindx + pos(i,1);
    yindxn = yindx + pos(i,2);
    xkeep = find(xindxn > 0 & xindxn <= sz(2));
    ykeep = find(yindxn > 0 & yindxn <= sz(1));
    discrep(i) = length(find(...
      mov(i).data(yindx(ykeep),xindx(xkeep)) ...
      - mov(1).data(yindxn(ykeep),xindxn(xkeep))));
  end
  