function disp(memm)
  s = struct(memm);
  protected_fields = {'chan2ind','offset','mm','d2v_slope','d2v_offset'};
  s = rmfield(s,protected_fields);
  s.channels = memm.header.channels;
  disp(s)
  