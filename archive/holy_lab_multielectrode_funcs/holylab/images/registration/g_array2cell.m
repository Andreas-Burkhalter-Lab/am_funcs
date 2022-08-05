function gc = g_array2cell(g)
  szg = size(g);
  n_dims = szg(end);
  colons = repmat({':'},1,n_dims);
  colons{end+1} = 1;
  gc = cell(1,n_dims);
  for i = 1:n_dims
    colons{end} = i;
    gc{i} = g(colons{:});
  end
end

