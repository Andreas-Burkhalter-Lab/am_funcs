//#define IMITERATORS_MAIN
#ifdef IMITERATORS_MAIN
#include <iostream>
#endif

// Iterator types:
//
//   pixIterator: a class to loop over all elements (pixels) in an
//     array by integer indexing.  It features edge-detection and
//     coordinate-extraction, as well as convenient methods for
//     defining sub-ranges of the array (useful for multithreading).
//     It is also written in a way that allows a single pixI to be
//     used for any number of arrays/images of the same size, so that you
//     can easily extract corresponding elements.  This is probably
//     the most broadly useful of all the classes.
//
//   colIterator: for accessing a 1-dimensional "pencil" ("column") of an image
//   multiDim1dIterator: for sequentially accessing all 1-dimensional
//     pencils oriented in a particular direction.
//   
//     These previous two classes are useful, for example, if you want
//     to filter an array along a particular dimension. (And if you
//     have a separable filter, you can apply the complete filter by
//     orienting the pencils along each dimension sequentially.)
//
//   nnIterator: nearest-neighbor iterator, for sequential access to all pixels
//     with nearest neighbors on all sides. You can effectively do the
//     same thing with pixIterator, and you should use pixIterator if
//     you need to do this while keeping two images "in sync" (this
//     was written first, which is why it's still around).



// This class provides access to a "vector" of values, in which the
// storage in memory may not be sequential but instead skip over a
// fixed number of elements. This situation arises when you want to
// look at a row or column of a matrix, for example.
template <class T>
class colIterator {
private:
  T* firstelement;
  T* endelement;
  T* curelement;
  size_t skipsize;
public:
  // Constructors
  colIterator() {
    firstelement = NULL;
    curelement = NULL;
    endelement = NULL;
    skipsize = 1;
  }
  colIterator(T* fe,T* ee,int ss) {
    firstelement = fe;
    curelement = fe;
    endelement = ee;
    skipsize = ss;
  }
  colIterator(const colIterator& ci) {
    firstelement = ci.firstelement;
    endelement = ci.endelement;
    curelement = ci.firstelement;
    skipsize = ci.skipsize;
  }
  
  // Pointer access operations
  T* begin() {return firstelement;}
  T* end() {return endelement;}
  bool operator<(T* comp) {return curelement < comp;}
  bool operator<=(T* comp) {return curelement <= comp;}
  bool operator>(T* comp) {return curelement > comp;}
  bool operator>=(T* comp) {return curelement >= comp;}

  // Value access operations
  T& operator*() {return *curelement;}
  T operator[](int offset) {return curelement[offset*skipsize];}

  // Increment/decrement operations
  colIterator& operator++(int) {curelement += skipsize; return *this;}
  colIterator& operator--(int) {curelement -= skipsize; return *this;}
  colIterator operator-(int);

  // Setting operations
  colIterator& operator=(T* ce) {curelement = ce; return *this;}

  // Information
  int length() {return (endelement-firstelement)/skipsize;}
};

template <class T>
colIterator<T> colIterator<T>::operator-(int offset)
{
  colIterator<T> ret(firstelement,endelement,skipsize);
  ret.curelement = curelement - offset*skipsize;
  return ret;
}

// This class provides a convenient way to iterate over all the rows
// or columns of a matrix, or more generally all 1-d "pencils" of a higher
// dimensional object.  It also supports ranging, which is useful for
// multithreaded applications in which you want to assign different
// chunks of the array to different processors---see the function
// mdi_split below.
// Constructor syntax:
//   multiDim1dIterator(T* first_element_pointer,int *dimension_size_vector,int number_of_dimensions,int chosen_dimension_index,int *rangestart,int *rangeend)
// See below (at end of this file) for usage examples
template <class T>
class multiDim1dIterator {
private:
  T* firstelement;       // Pointer to first element of array
  int n_dims;            // Number of dimensions in multidimensional array
  int *dims;             // Dimension sizes of multidimensional array
  int *range_start;      // "Upper left" corner of active region
  int *range_end;        // "Lower right" corner of active region
  int *current;          // Current coordinates (positioned at start of "column")
  size_t *cumSkip;       // Number of elements per increment in dimension
  int chosenDim;          // The chosen "column" dimension for 1d iterator
public:
  // Constructors & destructors
  multiDim1dIterator(T* fe,const int *d,int nd,int dI,const int *rs = NULL,const int *re = NULL);
  ~multiDim1dIterator();

  // Comparison
  bool at_end() {return current[n_dims] > 0;}  // uses beyond-the-end on dimensions

  // Increment
  multiDim1dIterator& inc_col();

  // Column extraction
  colIterator<T> col_iterator();
};

// A utility function for setting rangestart and rangeend given that
// you are going to split the job into a number of equal-sized
// chunks. Output is 0 on success and -1 on failure (which occurs for
// one-dimensional images)
// d = vector of dimensions sizes
// nd = # of dimensions
// dI = working dimension index (the pencil dimension)
// n_splits is the number of chunks desired
// split_index is the index [0 ... n_splits-1] of the current chunk
// rs and re are the output rangestart and rangeend
int mdi_split(const int *d,int nd,int dI,int n_splits,int split_index,int *rs,int *re)
{
  int dimIterator;
  for (dimIterator = 0; dimIterator < nd; dimIterator++) {
    rs[dimIterator] = 0;
    re[dimIterator] = d[dimIterator];
  }
  dimIterator = 0;        // split along the first dimension...
  if (dI == 0)
    dimIterator = 1;      // ...unless this is the dimension of the pencil
  if (dimIterator >= nd || split_index >= n_splits || split_index < 0)
    return -1;

  rs[dimIterator] = (d[dimIterator]*split_index)/n_splits;
  re[dimIterator] = (d[dimIterator]*(split_index+1))/n_splits;
  return 0;
}

template <class T>
multiDim1dIterator<T>::multiDim1dIterator(T* fe,const int *d,int nd,int dI,const int *rs,const int *re) {
  int dimIterator;

  firstelement = fe;
  chosenDim = dI;
  n_dims = nd;
  dims = new int [nd];
  range_start = new int [nd+1];      //+1 to handle incrementing a 1d mdi
  range_end = new int [nd+1];
  current = new int [nd+1];
  cumSkip = new size_t [nd+1];       // Last one contains total # elements

  for (dimIterator = 0; dimIterator < n_dims; dimIterator++) {
    dims[dimIterator] = d[dimIterator];
    if (rs != NULL)
      range_start[dimIterator] = rs[dimIterator];
    else
      range_start[dimIterator] = 0;   // default start at beginning
    if (re != NULL)
      range_end[dimIterator] = re[dimIterator];
    else
      range_end[dimIterator] = dims[dimIterator];   // default end at end
    current[dimIterator] = range_start[dimIterator];   // initialize to start of range
  }
  // Add a last ficticious dimension that gets incremented when we've
  // exhausted the range.  This keeps us from having to handle the 1d
  // case differently.
  range_start[n_dims] = 0;
  current[n_dims] = 0;
  range_end[n_dims] = 2;   // 2 will allow the fake dimension to be incremented

  *cumSkip = 1;
  for (dimIterator = 1; dimIterator <= n_dims; dimIterator++) {
    cumSkip[dimIterator] = cumSkip[dimIterator-1]*dims[dimIterator-1];
  }
}

template <class T>
multiDim1dIterator<T>::~multiDim1dIterator()
{
  delete [] dims;
  delete [] range_start;
  delete [] range_end;
  delete [] current;
  delete [] cumSkip;
}

template <class T>
multiDim1dIterator<T>& multiDim1dIterator<T>::inc_col()
{
  int dimIterator = 0;
  while (dimIterator <= n_dims) {
    if (dimIterator == chosenDim)
      dimIterator++;   // skip over the chosen column dimension
    else if (++current[dimIterator] >= range_end[dimIterator]) {
      current[dimIterator] = range_start[dimIterator];
      dimIterator++;
    }
    else
      break;  // we don't need to "carry" anymore
  }
}

template <class T>
colIterator<T> multiDim1dIterator<T>::col_iterator()
{
  int offset,len,dimSkip,dimIterator;

  // Calculate offset of first element in column
  offset = 0;
  for (dimIterator = 0; dimIterator < n_dims; dimIterator++)
    offset += cumSkip[dimIterator]*current[dimIterator];
  
  // Calculate the remaining items
  dimSkip = cumSkip[chosenDim];
  len = range_end[chosenDim]-range_start[chosenDim];

  // Return the column iterator
  return colIterator<T>(firstelement+offset,firstelement+offset+len*dimSkip,dimSkip);
}

// This class provides sequential access to all pixels with nearest
// neighbors on all sides; it also provides convenient access to those
// neighbors.  Note that pixIterator (below) can provide similar
// functionality; pixIterator is preferred when you want to do
// iteration over pixels in several images and keep them in sync
// (e.g., in image registration).
// To get the value of the current pixel:
//   *nni or nni[0]
// The pixel "to the left" (or, for matlab, "above"), i.e. with
// first coordinate one lower:
//   nni[-nni.dimSkip(0)]
// The pixel "below" (matlab: "to the right"), i.e., with second coordinate
// one higher:
//   nni[nni.dimSkip(1)]
// To move on to the next valid pixel:
//   nni++
// Constructor:
//   nnIterator(T* image_data_pointer,int *dimension_sizes,int number_of_dimensions,T* rangestart,T* rangeend)
// The last two are optional.  To divide an image into 4 pieces (say,
// if you're using 4CPUs), the third piece would be specified by
//    rangestart = image_data_pointer + npixels/2
//    rangeend = image_data_pointer + (3*npixels)/4
// Finally, note: dimensions with size 1 will be ignored.  This could
// be a possible source of confusion, since the object as indexed by an
// nnIterator may have smaller dimensionality than the nominal
// dimensionality passed to the constructor.  You can call the nDims()
// function to learn about how many dimensions were employed.
// Also note: if the array has dimension size 2 in any coordinate,
// then there are _no_ pixels that have nearest neighbors on all sides
// (all pixels are edge pixels), so this iterator will not give you
// any pixels.

template <class T>
class nnIterator {
private:
  T *firstelement;         // Pointer to first element of array
  int n_dims;              // Number of dimensions in multidimensional array
  int *dims;               // Dimension sizes of multidimensional array
  T *current;              // Pointer to current element
  int *current_c;          // Current coordinates
  const T *rangestart,*rangeend; // Region to work in
  size_t *cumSkip;         // Number of elements per increment in dimension
public:
  // Constructors & destructors
  nnIterator(T* fe,const int *d,int nd,const float *rs = NULL,const float *re = NULL);
  ~nnIterator();

  // Starting and ending
  void restart();    // Call this to position at first pixel in range
		     // (this is already called in the constructor)
  bool at_end() {return current >= rangeend;}

  // Increment
  nnIterator& operator++(int);

  // Access
  T operator[](int offset) {return current[offset];}
  T& operator*() {return *current;}

  // Dimensional information
  int nDims() {return n_dims;}
  size_t dimSkip(int dimensionIndex) {return cumSkip[dimensionIndex];}
};

template <class T>
nnIterator<T>::nnIterator(T* fe,const int *d,int nd,const float *rs,const float *re) {
  int dimIterator,dimIterator2;

  firstelement = fe;
  dims = new int [nd+1];
  current_c = new int [nd+1];        // +1 to handle final "carry"
  cumSkip = new size_t [nd+1];       // Last one contains total # elements

  // Copy over the dimension sizes, skipping dimensions with size 1
  for (dimIterator = 0,dimIterator2 = 0; dimIterator < nd; dimIterator++) {
    if (d[dimIterator] > 1) {
      dims[dimIterator2] = d[dimIterator];
      dimIterator2++;
    }
  }
  n_dims = dimIterator2;
  dims[n_dims] = 3;  // ficticious size of ficticious dimension

  *cumSkip = 1;
  for (dimIterator = 1; dimIterator <= n_dims; dimIterator++)
    cumSkip[dimIterator] = cumSkip[dimIterator-1]*dims[dimIterator-1];

  if (rs != NULL)
    rangestart = rs;
  else
    rangestart = firstelement;
  if (re != NULL)
    rangeend = re;
  else
    rangeend = firstelement+cumSkip[n_dims];

  // Set the initial coordinates & current pointer position
  restart();
}

template <class T>
nnIterator<T>::~nnIterator()
{
  delete [] dims;
  delete [] current_c;
  delete [] cumSkip;
}

template <class T>
void nnIterator<T>::restart()
{
  int elementIndex,dimIterator;
  // Set the initial coordinates; this may involve moving the pointer
  // ahead of rangestart to ensure that the first element has all of
  // its nearest neighbors
  elementIndex = rangestart - firstelement;
  for (dimIterator = n_dims-1; dimIterator >= 0; dimIterator--) {
    current_c[dimIterator] = elementIndex / cumSkip[dimIterator];
    elementIndex -= current_c[dimIterator] * cumSkip[dimIterator];
  }
  current_c[n_dims] = 0;   // ficticious beyond-the-end coordinate
  // Make the minimum in each coordinate 1, to ensure a left nearest-neighbor
  for (dimIterator = 0; dimIterator < n_dims; dimIterator++)
    if (current_c[dimIterator] < 1)
      current_c[dimIterator] = 1;
  // Make sure each coordinate is at least 1 away from the right
  // border, to ensure a right nearest-neighbor
  for (dimIterator = 0; dimIterator < n_dims; dimIterator++)
    if (current_c[dimIterator] >= dims[dimIterator]-1) {
      current_c[dimIterator] = 1;
      current_c[dimIterator+1]++;    // the carry operation
    }
	
  // Convert these coordinates to the current element pointer
  // Note the <= n_dims (rather than < n_dims); this is necessary to
  // ensure that we properly account for carrying into the
  // beyond-the-end coordinate.
  current = firstelement;
  for (dimIterator = 0; dimIterator <= n_dims; dimIterator++)
    current += current_c[dimIterator] * cumSkip[dimIterator];
}

template <class T>
nnIterator<T>& nnIterator<T>::operator++(int)
{
  int dimIterator;

  current_c[0]++;
  current++;
  if (current_c[0] >= dims[0]-1) {
    // We need to wrap around
    current_c[0] = 1;
    dimIterator = 1;
    while (dimIterator <= n_dims) {
      if (++current_c[dimIterator] >= dims[dimIterator]-1) {
	current_c[dimIterator] = 1;
	dimIterator++;
      } else
	break;      // all done carrying
    }
    // Recalculate the pointer position
    current = firstelement;
    for (dimIterator = 0; dimIterator <= n_dims; dimIterator++)
      current += current_c[dimIterator] * cumSkip[dimIterator];
  }
}

// pixIterator: a class to loop over all pixels in an image, using
// integer indexing.
//
// Construction:
//    pixI(d,nd) d is an int* array of dimension lengths, nd = # of dimensions
//    pixI(d,nd,false) use this if you want your coordinates to include
//      dimensions that were allocated as a singleton (default "true" is
//      to parametrize image as if "squeeze" had first been called on it)
//    pixI(d,nd,true,rs,re) allows you to specify starting and ending
//      indices that are a subset of the image, typically for use
//      in multithreading.  For example, if there are a total of n pixels
//      in the image, an iterator that would operate on the first of 4
//      chunks would have rs = 0, re = n/4; the second of 4 chunks would
//      have rs = n/4, re = n/2; and so forth.
//
// Usage: if im is a pointer to the image data (e.g., of type float*)
// and pixI is a pixIterator, then
//    im[pixI] gives the current pixel value
//    pixI.coord(2) is the coordinate of the 3rd dimension
//    pixI++ moves on to the next pixel
//    pixI.at_end() tests whether you've iterated beyond the region
//       covered by this iterator
//    pixI.restart() goes back to the beginning of the range
//    pixI.on_edge() is true if the current pixel is on the boundary of
//      the image
// There are other things you can do, see below.

class pixIterator {
private:
  long current;            // Memory offset of current pixel
#ifdef MAXDIMS
  int current_c[MAXDIMS+1];
  int dims[MAXDIMS+1];
  size_t cumSkip[MAXDIMS+1];
#else
  int *current_c;          // Coordinates of current pixel
  int *dims;               // Dimension sizes of multidimensional array
  size_t *cumSkip;         // Number of elements per increment in dimension
#endif
  int n_dims;              // Number of dimensions in multidimensional array
  long rangestart,rangeend;// Region to work in

  void allocate(int nd);
public:
  // Constructors, destructors, and related functions
  pixIterator(const int *d,int nd,bool skip_unity=true,long rs=0,long re=-1);
  pixIterator();  // empty constructor useful in arrays, etc.
  void initialize(const int *d,int nd,bool skip_unity=true,long rs=0,long re=-1);
  void setrange(long rs,long re);
  void copy(const pixIterator &pI);
  pixIterator(const pixIterator &pI) {copy(pI);}
  pixIterator& operator=(const pixIterator &pI) {copy(pI);}
  void clear();
  ~pixIterator();

  // Starting and ending
  void restart();    // Call this to position at first pixel in range
		     // (this is called in the constructor)
  bool at_end() const {return current >= rangeend;}

  // Increment current pixel
  pixIterator& operator++(int);
  pixIterator& inc_nocarry(int i=1);// increment without "wrapping" coordinates
  pixIterator& carry();             // to clean up after inc_nocarry

  // Current pixel info
  operator long() const {return current;}
  int coord(int dimIndex) const {return current_c[dimIndex];}
  bool on_edge() const;
  bool on_right_edge() const; // true if at "upper" edge
  bool on_edge(const bool *mask) const;  // set to true for the dimensions you care about
  bool on_right_edge(const bool *mask) const; // for specific dimensions
  bool on_edge_first() const;   // true if on an edge of the first coordinate

  // Dimensional information
  int nDims() const {return n_dims;}
  int dimSize(int dimIndex) const {return dims[dimIndex];}
  size_t dimSkip(int dimensionIndex) const {return cumSkip[dimensionIndex];}
  size_t numel() const {return cumSkip[n_dims];}

  // Offset information
  size_t offset_from_coords(const pixIterator &pI) const;
} __attribute__((aligned(64)));  // for use in multithreading to ensure diff cache lines (this makes a big difference)

void pixIterator::allocate(int nd)
{
#ifndef MAXDIMS
  dims = new int [nd+1];             // +1 to handle final "carry" (sentinel)
  current_c = new int [nd+1];        // +1 to handle final "carry" (sentinel)
  cumSkip = new size_t [nd+1];       // Last one contains total # elements
#endif
}  

pixIterator::pixIterator(const int *d,int nd,bool skip_unity,long rs,long re)
{
#ifndef MAXDIMS
  dims = NULL;
  current_c = NULL;
  cumSkip = NULL;
#endif
  initialize(d,nd,skip_unity,rs,re);
}

pixIterator::pixIterator()
{
#ifndef MAXDIMS
  dims = NULL;
  current_c = NULL;
  cumSkip = NULL;
#endif
  n_dims = 0;
  current = 0;
  rangestart = 0;
  rangeend = -1;
}

void pixIterator::initialize(const int *d,int nd,bool skip_unity,long rs,long re)
{
  int dimIterator,dimIterator2;

  clear();
  allocate(nd);

  // Copy over the dimension sizes, optionally skipping dimensions with size 1
  for (dimIterator = 0,dimIterator2 = 0; dimIterator < nd; dimIterator++) {
    if (!skip_unity || d[dimIterator] > 1) {
      dims[dimIterator2] = d[dimIterator];
      dimIterator2++;
    }
  }
  n_dims = dimIterator2;
  dims[n_dims] = 3;  // ficticious size of ficticious dimension (for sentinel)

  *cumSkip = 1;
  for (dimIterator = 1; dimIterator <= n_dims; dimIterator++)
    cumSkip[dimIterator] = cumSkip[dimIterator-1]*dims[dimIterator-1];

  setrange(rs,re);
}

void pixIterator::setrange(long rs,long re)
{
  rangestart = rs;
  if (re != -1)
    rangeend = re;
  else
    rangeend = cumSkip[n_dims];
  // Set the initial coordinates & current pointer position
  restart();
}

void pixIterator::copy(const pixIterator &pI)
{
  clear();
  
  current = pI.current;
  n_dims = pI.n_dims;
  rangestart = pI.rangestart;
  rangeend = pI.rangeend;
  allocate(n_dims);

  int dimIterator;
  for (dimIterator = 0; dimIterator <= n_dims; dimIterator++) {
    current_c[dimIterator] = pI.current_c[dimIterator];
    dims[dimIterator] = pI.dims[dimIterator];
    cumSkip[dimIterator] = pI.cumSkip[dimIterator];
  }
}
    

pixIterator::~pixIterator()
{
  clear();
}

void pixIterator::clear()
{
#ifndef MAXDIMS
  if (dims != NULL)
    delete [] dims;
  if (current_c != NULL)
    delete [] current_c;
  if (cumSkip != NULL)
    delete [] cumSkip;
  dims = NULL;
  current_c = NULL;
  cumSkip = NULL;
#endif

  n_dims = 0;
  current = 0;
  rangestart = 0;
  rangeend = -1;
}

void pixIterator::restart()
{
  int dimIterator,remainder;

  current = rangestart;
  remainder = current;
  for (dimIterator = n_dims; dimIterator >= 0; dimIterator--) {
    current_c[dimIterator] = remainder / cumSkip[dimIterator];
    remainder -= current_c[dimIterator] * cumSkip[dimIterator];
  }
}

pixIterator& pixIterator::operator++(int)
{
  int dimIterator;

  current++;
  // Increment coordinates, with carry. Note that this carry will be
  // complete only if we're up-to-date, i.e., any inc_nocarry calls
  // have been cleaned up by calling carry.
  dimIterator = 0;
  while (dimIterator <= n_dims) {
    if (++current_c[dimIterator] >= dims[dimIterator]) {
      current_c[dimIterator] -= dims[dimIterator];
      dimIterator++;
    } else
      break;      // all done carrying
  }
}

pixIterator& pixIterator::inc_nocarry(int i)
{
  int dimIterator;

  current += i;
  // Increment coordinates, without carry. Make sure you eventually call carry!
  current_c[0] += i;
}

pixIterator& pixIterator::carry()
{
  // This is a "complete" carry, i.e., it works no matter how far off you are
  int dimIterator;

  dimIterator = 0;
  while (dimIterator < n_dims) {
    if (current_c[dimIterator] >= dims[dimIterator]) {
      current_c[dimIterator+1] += current_c[dimIterator]/dims[dimIterator];
      current_c[dimIterator] = current_c[dimIterator] % dims[dimIterator];
      dimIterator++;
    } else
      break;      // all done carrying
  }
}

// Calculate the memory offset from the coordinates of another pixIterator
size_t pixIterator::offset_from_coords(const pixIterator &pI) const
{
  size_t ret;
  int dimIterator;

  for (ret = 0, dimIterator = 0; dimIterator < n_dims; dimIterator++)
    ret += pI.coord(dimIterator) * cumSkip[dimIterator];

  return ret;
}

// The following inlines have been tested and found to improve execution speed
inline bool pixIterator::on_edge() const
{
  const int *c,*d;
  const int *dend = dims + n_dims;
  
  for (c = current_c, d = dims; d != dend; c++,d++)
    if (*c == 0 || *c == *d-1)
      return true;

  return false;
}

inline bool pixIterator::on_right_edge() const
{
  const int *c,*d;
  const int *dend = dims + n_dims;
  
  for (c = current_c, d = dims; d != dend; c++,d++)
    if (*c == *d-1)
      return true;

  return false;
}

inline bool pixIterator::on_edge(const bool *mask) const
{
  const int *c,*d;
  const int *dend = dims + n_dims;
  
  for (c = current_c, d = dims; d != dend; c++,d++,mask++)
    if (*mask && (*c == 0 || *c == *d-1))
	return true;

  return false;
}

inline bool pixIterator::on_right_edge(const bool *mask) const
{
  const int *c,*d;
  const int *dend = dims + n_dims;
  
  for (c = current_c, d = dims; d != dend; c++,d++,mask++)
    if (*mask && *c == *d-1)
	return true;

  return false;
}

inline bool pixIterator::on_edge_first() const
{
  return (current_c[0] == 0 || current_c[0] == dims[0]-1);
}




// Functions below are for testing purposes, and also give an example
// of usage
#ifdef IMITERATORS_MAIN
void dimfunc(float x[],int *dims,int n_dims)
{
  int dimIndex,dimIndex2,n_elements;
  int *rangestart,*rangeend;

  rangestart = new int [n_dims];
  rangeend = new int [n_dims];

  n_elements = dims[0];
  for (dimIndex = 1; dimIndex < n_dims; dimIndex++)
    n_elements *= dims[dimIndex];

  for (dimIndex = 0; dimIndex < n_dims; dimIndex++) {
    // Full multidimensional iterator
    multiDim1dIterator<float> mdi(x,dims,n_dims,dimIndex);
    std::cout << "mdIterator in dimension " << dimIndex << ":\n";
    for (; !mdi.at_end(); mdi.inc_col()) {
      colIterator<float> ci = mdi.col_iterator();
      for (; ci < ci.end(); ci++)
	std::cout << *ci << ' ';
      std::cout << std::endl;
    }

    // First-half mditerator
    if (mdi_split(dims,n_dims,dimIndex,2,0,rangestart,rangeend) == 0) {
      multiDim1dIterator<float> mdif(x,dims,n_dims,dimIndex,rangestart,rangeend);
      std::cout << "First-half mdIterator in dimension " << dimIndex << ":\n";
      for (; !mdif.at_end(); mdif.inc_col()) {
	colIterator<float> ci = mdif.col_iterator();
	for (; ci < ci.end(); ci++)
	  std::cout << *ci << ' ';
	std::cout << std::endl;
      }
    }

    // Second-half mditerator
    if (mdi_split(dims,n_dims,dimIndex,2,1,rangestart,rangeend) == 0) {
      multiDim1dIterator<float> mdis(x,dims,n_dims,dimIndex,rangestart,rangeend);
      std::cout << "Second-half mdIterator in dimension " << dimIndex << ":\n";
      for (; !mdis.at_end(); mdis.inc_col()) {
	colIterator<float> ci = mdis.col_iterator();
	for (; ci < ci.end(); ci++)
	  std::cout << *ci << ' ';
	std::cout << std::endl;
      }
    }
  }

  std::cout << "Nearest-neighbor iterator:\n";
  nnIterator<float> nni(x,dims,n_dims);
  while (!nni.at_end()) {
    std::cout << *nni << ' ';
    nni++;
  }
  std::cout << std::endl;  

  std::cout << "First-half nearest-neighbor iterator:\n";
  nnIterator<float> nnif(x,dims,n_dims,x,x+n_elements/2);
  while (!nnif.at_end()) {
    std::cout << *nnif << ' ';
    nnif++;
  }
  std::cout << std::endl;  

  std::cout << "Second-half nearest-neighbor iterator:\n";
  nnIterator<float> nnis(x,dims,n_dims,x+n_elements/2);
  while (!nnis.at_end()) {
    std::cout << *nnis << ' ';
    nnis++;
  }
  std::cout << std::endl;  

  delete [] rangestart;
  delete [] rangeend;
}

void pixPrint(pixIterator& pixI)
{
  int dimIterator;

  for (dimIterator = 0; dimIterator < pixI.nDims(); dimIterator++)
    std::cout << pixI.coord(dimIterator) << ' ';
  std::cout << " (" << pixI << ")  on edge:" << pixI.on_edge() << std::endl;
}

int main()
{
  int skip;
  float x[] = {1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1,2.2};
  int len = 12;
  int dims[3];
  int n_dims,dimIndex;
  

  colIterator<float> ci1(x,x+len,1);
  colIterator<float> ci2(x,x+len,2);

  for (; ci1 < ci1.end(); ci1++)
    std::cout << *ci1 << ' ';
  std::cout << std::endl;

  for (ci2 = ci2.end()-1; ci2 >= ci2.begin(); ci2--)
    std::cout << *ci2 << ' ';
  std::cout << std::endl;

  ci2 = ci2.begin();
  std::cout << ci2[2] << ' ' << ci1[-5] << std::endl;

  for (skip = 1; skip < 4; skip++) {
    colIterator<float> ci(x,x+len,skip);
    for (; ci < ci.end(); ci++)
      std::cout << *ci << ' ';
    std::cout << std::endl;
  }

  std::cout << "Representing as a 1-d object" << std::endl;
  dims[0] = len;
  dimfunc(x,dims,1);
      
  std::cout << "Representing as a 2-d object" << std::endl;
  dims[0] = 4;
  dims[1] = 3;
  dimfunc(x,dims,2);

  std::cout << "Represent as a 3-d object" << std::endl;
  dims[0] = 2;
  dims[1] = 2;
  dims[2] = 3;
  dimfunc(x,dims,3);

  std::cout << "\n\n\npixIterator exercises" << std::endl;
  dims[0] = 4;
  dims[1] = 3;
  //pixIterator pixI(dims,2,5,9);
  pixIterator pixI(dims,2);
  for (; !pixI.at_end(); pixI++)
    pixPrint(pixI);

  return 0;
}
#endif
