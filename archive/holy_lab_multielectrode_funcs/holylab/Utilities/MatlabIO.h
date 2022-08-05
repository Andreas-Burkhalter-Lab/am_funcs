#ifndef MatlabIO_h
#define MatlabIO_h

#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>
#include <stdint.h>
#include "MatlabTraits.h"

/** \brief Support for converting Matlab matrices to/from C++ types and reading/writing .mat files
 * \details You should include "mex.h" (for a MEX file) and/or "mat.h" (if you want to read/write .mat files) before including this file.  Use mat.h if your goal is to create a stand-alone program that can read/write Matlab data.\n
 * Likewise, for Eigen support be sure to include Eigen/Core before including this file.\n
 * Stand-alone applications should be compiled as
 * \verbatim mex -f $MATLABINSTALLDIR/bin/matopts.sh cppfiles -o executable \endverbatim
 * You may encounter problems of several different types, most of which can be fixed by editing the matops.sh file.  Make sure you make the modifications in the section appropriate to your architecture.
 * - Problems with //-style comments: delete the -ansi flag from CFLAGS
 * - Problems linking to the correct libraries: find the variable RPATH and change the part that says
 * \verbatim -rpath-link,$TMW_ROOT/bin/$Arch \endverbatim
 * to \verbatim -rpath=$TMW_ROOT/bin/$Arch \endverbatim
 * Recompile and try again.
 */
namespace MatlabIO {
//
// std::vector wrappers and support functions
//
/** \brief Convert an mxArray (Matlab) vector into a std::vector
 *  \details Because std::vectors can be resized---and this would wreak havoc with the Matlab memory manager---this makes a copy of the data in the mxArray.  Use the Eigen types if you want to avoid copies.
 *  \param pMxA A pointer to the matlab mxArray
 *  \param x A std::vector variable (e.g., std::vector<double>)
 *  \exception std::runtime_error is thrown if the mxArray's type does not match that of x, or if the mxArray has more than 1 dimension.
 */
template <typename T>
void mx2vector(const mxArray* pMxA,std::vector<T> &x)
{
  // Check that the mxArray is of the right type
  if (!MatlabTraits<T>::IsA(pMxA))
    throw std::runtime_error(std::string("mx2vector: mxArray is not of type ") +
			     std::string(MatlabTraits<T>::ComplexName()) +
			     std::string(" ") +
			     std::string(MatlabTraits<T>::RealName()));
  // Check that the mxArray is a vector. We don't insist that it be
  // row/column/whatever, only that it have only one dimension with
  // nontrivial "length"
  int n_dims = mxGetNumberOfDimensions(pMxA);
  const mwSize *sz = mxGetDimensions(pMxA);
  int n_bigger_than_one = 0;
  for (int i = 0; i < n_dims; i++)
    if (sz[i] > 1)
      n_bigger_than_one++;
  if (n_bigger_than_one > 1)
    throw std::runtime_error("mx2vector: mxArray must be a vector");
  // Perform the conversion
  x.erase(x.begin(),x.end());
  T* pData = (T*) mxGetData(pMxA);
  x.insert(x.begin(),pData,pData+mxGetNumberOfElements(pMxA));
}


/** \brief Convert a std::vector into an mxArray (Matlab) row vector
 *  \param x A std::vector
 *  \return A pointer to the created mxArray
 *  \note Because Matlab needs to manage its own memory, this makes a copy of the data in x.
 */
template <typename T>
mxArray* vector2mx(const std::vector<T> &x)
{
  int sz[2];
  mxArray *pMxA;
  T* pData;

  sz[0] = 1;
  sz[1] = x.size();
  pMxA = MatlabTraits<T>::mxAllocator(2,sz);
  pData = (T*) mxGetData(pMxA);
  std::copy(x.begin(),x.end(),pData);
  return pMxA;
}



#ifdef EIGEN_CORE_H
//
//  Eigen wrappers and conversion functions
//
template <typename T>
struct MapMx {
  typedef Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> MatrixType;
  typedef Eigen::Matrix<T,Eigen::Dynamic,1>              ColVectorType;
  typedef Eigen::Map<MatrixType>                         MapMatrixType;
  typedef Eigen::Map<const MatrixType>                   MapConstMatrixType;
  typedef Eigen::Map<ColVectorType>                      MapColVectorType;
  typedef Eigen::Map<const ColVectorType>                MapConstColVectorType;
};
/** \brief Double matrix type
 */
typedef MapMx<double>::MapMatrixType MapMxD;
/** \brief Float matrix type
 */
typedef MapMx<float>::MapMatrixType MapMxF;
/** \brief Int8 matrix type
 */
typedef MapMx<int8_t>::MapMatrixType MapMxInt8;
/** \brief Uint8 matrix type
 */
typedef MapMx<uint8_t>::MapMatrixType MapMxUint8;
/** \brief Int16 matrix type
 */
typedef MapMx<int16_t>::MapMatrixType MapMxInt16;
/** \brief Uint16 matrix type
 */
typedef MapMx<uint16_t>::MapMatrixType MapMxUint16;
/** \brief Int32 matrix type
 */
typedef MapMx<int32_t>::MapMatrixType MapMxInt32;
/** \brief Uint32 matrix type
 */
typedef MapMx<uint32_t>::MapMatrixType MapMxUint32;
/** \brief Int64 matrix type
 */
typedef MapMx<int64_t>::MapMatrixType MapMxInt64;
/** \brief Uint64 matrix type
 */
typedef MapMx<uint64_t>::MapMatrixType MapMxUint64;

/** \brief Convert an mxArray (Matlab matrix) into appropriate Eigen type
 *  \details This doesn't copy the data, this just "wraps" the Matlab-allocated memory as an Eigen array
 *  \param pMxA A pointer to the matlab mxArray
 *  \param X An Eigen::Map variable (e.g., MapMxD)
 *  \exception std::runtime_error is thrown if the mxArray's type does not match that of X, or if the mxArray has more than 2 dimensions.
 */
template <typename Derived>
void mx2eigen(const mxArray* pMxA,Eigen::Map<Derived> &X)
{
  typedef typename Derived::Scalar Scalar;
  // Check that the mxArray is of the right type
  if (!MatlabTraits<Scalar>::IsA(pMxA))
    throw std::runtime_error(std::string("mx2eigen: mxArray is not of type ") +
			     std::string(MatlabTraits<Scalar>::ComplexName()) +
			     std::string(" ") +
			     std::string(MatlabTraits<Scalar>::RealName()));
  if (mxGetNumberOfDimensions(pMxA) > 2)
    throw std::runtime_error("mx2eigen: mxArray cannot have more than 2 dimensions");
  // Perform the conversion
  new (&X) Eigen::Map<Derived>((Scalar*) mxGetData(pMxA), mxGetM(pMxA), mxGetN(pMxA));
}

/** \brief Convert an Eigen type into an mxArray (Matlab matrix)
 *  \param X An Eigen type (Matrix, Array, Map, ...)
 *  \return A pointer to the created Matlab mxArray
 *  \note Because Matlab needs to manage its own memory, this makes a copy of the data in X.  As a zero-copy alternative, you can first allocate the mxArray, and then use mx2eigen to get an Eigen type that maps the allocated memory.  Example:
 *  \code
 *  // Create the output array
 *  int sz[] = {3,5};
 *  mxArray *outputMx = MatlabTraits<double>::mxAllocator(2,sz); // alternatively use mxCreateDoubleMatrix
 *  // Map to an Eigen matrix
 *  MapMxD outputEigen(NULL,0,0);
 *  mx2eigen(outputMx,outputEigen);
 *  // Now "fill" the memory with the result of your calculation
 *  myFunction(input1,...,outputEigen)
 *  \endcode
 */
template <typename Derived>
mxArray* eigen2mx(const Eigen::DenseBase<Derived> &X)
{
  typedef typename Derived::Scalar Scalar;
  typedef typename MapMx<Scalar>::MapMatrixType MapMxType;
  int sz[2];
  mxArray *pMxA;

  sz[0] = X.rows();
  sz[1] = X.cols();
  pMxA = MatlabTraits<Scalar>::mxAllocator(2,sz);
  // Because X might be an expression type, we need to evaluate the
  // coefficients explicitly rather than using std::copy
  MapMxType Xtmp(NULL,0,0);
  mx2eigen(pMxA,Xtmp);
  Xtmp = X.eval();
  /*
  Scalar *pData = (Scalar *) mxGetData(pMxA);
  for (int i = 0; i < X.size(); i++)
    pData[i] = X(i);
  */
  return pMxA;
}
#endif   // EIGEN_CORE_H

// For the next functions we have to make certain that the .mat file
// I/O is defined
#if defined(mat_h) || defined(mat_published_c_api_h)

/** \class Load
 *  \brief Load data from .mat files
 */
class Load {
public:
// LIFECYCLE
  /** Constructor
   *  \param filename A char* containing the name of the file to open
   *  \exception std::runtime_error is thrown if the file cannot be opened
   */
  Load(const char* filename) : mMatFP(matOpen(filename,"r")) {
    if (mMatFP == NULL)
      throw std::runtime_error("file open failure");
    // Read the names of all of the variables in the file, to support queries
    mxArray* pMxTmp;
    const char* name;
    pMxTmp = matGetNextVariableInfo(mMatFP,&name);
    while (pMxTmp != NULL) {
      mVariableList.push_back(name);
      pMxTmp = matGetNextVariableInfo(mMatFP,&name);
    }
  }
  /** Destructor (closes the file)
   *  \exception std::runtime_error is thrown if the file cannot be closed
   */
  ~Load() {
    if (matClose(mMatFP) != 0)
      throw std::runtime_error("failure to close file");
  }
  // Copy constructor and assignment declared private to prevent their use

// OPERATIONS    
  /** Load named Matlab variable and return in Matlab "mx" format.
   *  This is here to handle cases where the data cannot be imported as an Eigen type (e.g., is not a matrix).
   *  \param varname A char* containing the name of a variable in the file
   *  \return A pointer to an an mxArray
   *  \exception std::runtime_error is thrown if the named variable cannot be read
   */
  mxArray* load(const char* varname) {
    read(varname);
    return mpMxA;
  }
#ifdef EIGEN_CORE_H
  /** Load named Matlab variable and convert to Eigen type
   *  \param varname A char* containing the name of a variable in the file
   *  \param X An Eigen::Map type (e.g., MapMxD)
   *  \exception std::runtime_error is thrown if the named variable cannot be read
   */
  template <typename Derived> void load(const char* varname,Eigen::Map<Derived> &X) {
    typedef typename Derived::Scalar Scalar;
    read(varname);
    // Since we have the name of the variable, we can use it in
    // error-reporting, so do our own type checking rather than using
    // mx2eigen
    if (!MatlabTraits<Scalar>::IsA(mpMxA))
      throw std::runtime_error(std::string(varname) +
			       std::string(" is not a variable of type ") +
			       std::string(MatlabTraits<Scalar>::ComplexName()) +
			       std::string(" ") +
			       std::string(MatlabTraits<Scalar>::RealName()));
    if (mxGetNumberOfDimensions(mpMxA) > 2)
      throw std::runtime_error(std::string(varname) + std::string(" has more than 2 dimensions"));
    new (&X) Eigen::Map<Derived>((Scalar*) mxGetData(mpMxA), mxGetM(mpMxA), mxGetN(mpMxA));
  }
#endif   // EIGEN_CORE_H

  /** Load a C++ scalar as a named variable
   *  \param s A scalar value (e.g., int, double, etc)
   *  \param varname A char* containing the name of the variable to load
   *  \exception std::runtime_error if write failure
   */
  template <typename T> void loadScalar(const char* varname,T& s) {
    read(varname);
    if (mxGetNumberOfElements(mpMxA) != 1)
      throw std::runtime_error(std::string(varname) + std::string(" is not a scalar"));
    s = T(mxGetScalar(mpMxA));
  }


  // INQUIRY
  /** Determine whether a particular named variable is in the file
   * \param name A string or character array with the variable name
   * \return true if the named variable is stored in the file
   */
  bool hasVariable(const char* name) {
    std::vector<std::string>::iterator vsi;
    vsi = find(mVariableList.begin(),mVariableList.end(),name);
    return (vsi != mVariableList.end());
  }


private:
  mxArray *mpMxA;
  MATFile *mMatFP;
  std::vector<std::string> mVariableList;
  
  // Read an mxArray into internal store
  void read(const char* varname) {
    mpMxA = matGetVariable(mMatFP,varname);
    if (mpMxA == NULL)
      throw std::runtime_error(std::string("Failure to read variable ") + std::string(varname));
  }
  // Declare these private to prevent copying
  Load(const Load &);
  Load & operator= (const Load &);
};



/** \brief Save data to .mat files
 */
class Save {
public:
// LIFECYCLE
  /** Constructor
   *  \param filename A char* containing the name of the file to open, or NULL if you do not actually want to open a file
   *  \exception std::runtime_error is thrown if the file cannot be opened
   */
  Save(const char* filename) {
    mMatFP = NULL;
    if (filename != NULL) {
      mMatFP = matOpen(filename,"w");
      if (mMatFP == NULL)
	throw std::runtime_error("file open failure");
    }
  }
  /** Destructor (closes the file)
   *  \exception std::runtime_error is thrown if the file cannot be closed
   */
  ~Save() {
    if (mMatFP != NULL) {
      if (matClose(mMatFP) != 0)
	throw std::runtime_error("failure to close file");
    }
  }
  // Copy constructor and assignment declared private to prevent their use

// OPERATIONS    
  /** Save Matlab array as named variable
   *  \param pMxA A pointer to an an mxArray
   *  \param varname A char* containing the name you want to use for the saved variable
   *  \exception std::runtime_error if write failure
   */
  void save(const mxArray* pMxA,const char* varname) {
    write(pMxA,varname);
  }
  /** Convert std::vector type to mxArray and save as named variable
   *  \param x A vector (std::vector)
   *  \param varname A char* containing the name you want to use for the saved variable
   *  \exception std::runtime_error if write failure
   */
  template <typename T> void save(const std::vector<T> &X,const char* varname) {
    mxArray *pMxA = vector2mx(X);
    write(pMxA,varname);
    mxDestroyArray(pMxA);
  }
#ifdef EIGEN_CORE_H
  /** Convert Eigen type to mxArray and save as named variable
   *  Since this involves making a temporary, a possibly-better approach is to first allocate the output as an mxArray, map to an Eigen type, and then fill the data in this type with the results of your calculation (see help to eigen2mx). Then save the mxArray.
   *  \param X An Eigen data type
   *  \param varname A char* containing the name you want to use for the saved variable
   *  \exception std::runtime_error if write failure
   */
  template <typename Derived> void save(const Eigen::DenseBase<Derived> &X,const char* varname) {
    mxArray *pMxA = eigen2mx(X);
    write(pMxA,varname);
    mxDestroyArray(pMxA);
  }
#endif  // EIGEN_CORE_H
  /** Save a C++ scalar as a named variable
   *  \param s A scalar value (e.g., int, double, etc)
   *  \param varname A char* containing the name you want to use for the saved variable
   *  \exception std::runtime_error if write failure
   */
  template <typename T> void saveScalar(T s,const char* varname) {
    const int sz[] = {1,1};
    mxArray *pMxA = MatlabTraits<T>::mxAllocator(2,sz);
    T* pData = (T*) mxGetData(pMxA);
    *pData = s;
    write(pMxA,varname);
  }

private:
  MATFile *mMatFP;
  
  void write(const mxArray* pMxA,const char* varname) {
    if (mMatFP == NULL)
      throw std::runtime_error("File is not open for writing");
    int status = matPutVariable(mMatFP,varname,pMxA);
    if (status != 0) {
      matClose(mMatFP);
      throw std::runtime_error(std::string("Failure to write variable ") + std::string(varname));
    }
  }
  // Declare these private to prevent copying
  Save(const Save &);
  Save & operator= (const Save &);
};


#endif // ifdef mat_h
}

#endif // MatlabIO_h
