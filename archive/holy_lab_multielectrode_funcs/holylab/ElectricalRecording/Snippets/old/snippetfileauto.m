function snippetfileauto(infilename,varargin)
  [fid,message] = fopen(infilename,'r','b');
  if (fid == -1)
    error(message);
  end
  h = ReadAIHeader(fid);
  scanrange = [0 h.nscans-1];
  channels = h.channels;
  chanIndx = find(channels > 2);
  channels = channels(chanIndx);
  fseek(fid,h.headersize,'bof');
  w = fread(fid,[h.numch 50000],'*int16');
  so = sfaparse(varargin);
  so.Fs = h.scanrate;
  so = snipoptions(so);
  wf = filterint16(so.condfiltb,so.condfilta,w,chanIndx);
  mna = mean(abs(double(wf)),2);
  thresh = 6*mna;
  if (length(varargin) > 0)
    snippetfile(fid,scanrange,channels,thresh,varargin{:});
  else
    snippetfile(fid,scanrange,channels,thresh);
  end


function so = sfaparse(P)
  iss = zeros(1,length(P));
  for i = 1:length(P)
    iss(i) = isstruct(P{i});
  end
  if (length(find(iss)) > 1)
    error('Only one structure input');
  end
  if any(iss)
    so = P{find(iss)};
  else
    so = [];
  end