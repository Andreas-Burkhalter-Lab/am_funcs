template <typename valType,typename intType>
void heapsort_sift_down(valType *val, intType *index, const intType l, const intType r)
{
  valType valtmp;
  intType indextmp, position, positionOld;

  valtmp = val[l];
  indextmp = index[l];
  positionOld = l;
  position = 2*l+1;
  while (position <= r) {
    if (position < r && val[position] < val[position+1])
      position++;
    if (valtmp >= val[position])
      break;
    val[positionOld] = val[position];
    index[positionOld] = index[position];
    positionOld = position;
    position = 2*position+1;
  }
  val[positionOld] = valtmp;
  index[positionOld] = indextmp;
}

template <typename valType,typename intType>
void heapsort_index(valType *val,intType *index,intType n)
{
  valType valtmp;
  intType indextmp, position;

  // Fill the index with 0:n-1
  for (position = 0; position < n; position++)
    index[position] = position;
  // Build the heap
  for (position = n/2-1; position >= 0; position--)
    heapsort_sift_down(val,index,position,n-1);
  // Order the elements
  for (position = n-1; position > 0; position--) {
    valtmp = val[0];
    indextmp = index[0];
    val[0] = val[position];
    index[0] = index[position];
    val[position] = valtmp;
    index[position] = indextmp;
    heapsort_sift_down(val,index,0,position-1);
  }
}
