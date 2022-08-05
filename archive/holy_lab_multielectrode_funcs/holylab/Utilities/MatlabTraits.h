/** \class MatlabTraits
 * \brief Information to facilitate working with different Matlab data types
 *
 * \details Matlab's external interface to C, with different functions for different data types, makes it somewhat awkward to use C++ templates.  This traits structure makes this easier.
 *
 * \note You'll need to \#include "mex.h" or \#include "mat.h" (depending on whether you are creating a MEX file or a stand-alone application) to use this file.  See MatEigenIO.h for notes on compiling such programs.
 * \param T The C++ type for which you want "access" to the corresponding Matlab support routines.  Complex numbers may be handled by, e.g., std::complex<float>.
 *
 * The most useful elements are the following:
 *   \li A function IsA(const mxArray *pA), which evaluates to true if the mxArray pointed to by pA is of the appropriate kind to store data of type T.
 *   \li A function mxAllocator(mwSize ndim,const mwSize *dims) which allocates an mxArray of the appropriate type to store data of the indicated size.
 *   \li RealName(), a function returning a character string naming the Real data type (useful for error-reporting).
 *   \li ComplexName(), a function returning either "real" or "complex" for numeric types as appropriate. For char this is "".
 *
 * The following are also available:
 *   \li A typedef Real, corresponding to the data type used to store values.  For most types this is just T; but if T is a \c std::complex<U>, this would be U.
 *   \li A constant IsNumeric, evaluates to 1 (true) if T is a numeric type.
 *   \li A constant IsComplex, evaluates to 1 (true) if T is complex.
 *   \li A function ClassID(), returning the Matlab enum value encoding the type.
 *   \li A function Complexity(), either mxREAL or mxCOMPLEX as appropriate to the complexity type.
 *   \li A function ScalarTypeCheck(const mxArray *pA), which checks only the scalar data type.  This differs from IsA in that IsA also checks the complexity flag.
 *
 * Example:
 * \include test_MatlabTraits.cpp
 * Output: \verbinclude test_MatlabTraits.out
 */

// Acknowledgments: this is inspired by Eigen::NumTraits

#ifndef MatlabTraits_h
#define MatlabTraits_h

#include <complex>
#include <stdint.h>

#ifndef DOXYGEN_SHOULD_SKIP_THIS
// Default implementation is for a real double
template<typename T>
struct GenericBaseMatlabTraits
{
  typedef T Real;
  enum {
    IsNumeric = 1,
    IsComplex = 0
  };
  // These next two are defined as functions rather than enums because
  // they need to have the appropriate type
  inline static mxClassID ClassID() { return mxDOUBLE_CLASS; }
  // The next is useful for error-reporting
  inline static const char* RealName() { return "double"; }
};

template<typename T>
struct BaseMatlabTraits : GenericBaseMatlabTraits<T> {};

// Specialize for other types by overriding various fields
// Version for float = single
template<>
struct BaseMatlabTraits<float> : GenericBaseMatlabTraits<float>
{
  inline static mxClassID ClassID() { return mxSINGLE_CLASS; }
  inline static const char* RealName() { return "single"; }
};

// Version for int8
template<>
struct BaseMatlabTraits<int8_t> : GenericBaseMatlabTraits<int8_t>
{
  inline static mxClassID ClassID() { return mxINT8_CLASS; }
  inline static const char* RealName() { return "int8"; }
};

// Version for uint8
template<>
struct BaseMatlabTraits<uint8_t> : GenericBaseMatlabTraits<uint8_t>
{
  inline static mxClassID ClassID() { return mxUINT8_CLASS; }
  inline static const char* RealName() { return "uint8"; }
};

// Version for int16
template<>
struct BaseMatlabTraits<int16_t> : GenericBaseMatlabTraits<int16_t>
{
  inline static mxClassID ClassID() { return mxINT16_CLASS; }
  inline static const char* RealName() { return "int16"; }
};

// Version for uint16
template<>
struct BaseMatlabTraits<uint16_t> : GenericBaseMatlabTraits<uint16_t>
{
  inline static mxClassID ClassID() { return mxUINT16_CLASS; }
  inline static const char* RealName() { return "uint16"; }
};

// Version for int32
template<>
struct BaseMatlabTraits<int32_t> : GenericBaseMatlabTraits<int32_t>
{
  inline static mxClassID ClassID() { return mxINT32_CLASS; }
  inline static const char* RealName() { return "int32"; }
};

// Version for uint32
template<>
struct BaseMatlabTraits<uint32_t> : GenericBaseMatlabTraits<uint32_t>
{
  inline static mxClassID ClassID() { return mxUINT32_CLASS; }
  inline static const char* RealName() { return "uint32"; }
};

// Version for int64
template<>
struct BaseMatlabTraits<int64_t> : GenericBaseMatlabTraits<int64_t>
{
  inline static mxClassID ClassID() { return mxINT64_CLASS; }
  inline static const char* RealName() { return "int64"; }
};

// Version for uint64
template<>
struct BaseMatlabTraits<uint64_t> : GenericBaseMatlabTraits<uint64_t>
{
  inline static mxClassID ClassID() { return mxUINT64_CLASS; }
  inline static const char* RealName() { return "uint64"; }
};

// Version for complex types
template<typename _Real>
struct BaseMatlabTraits<std::complex<_Real> > : GenericBaseMatlabTraits<_Real>
{
  typedef _Real Real;
  enum {
    IsComplex = 1
  };
  inline static mxClassID ClassID() { return BaseMatlabTraits<_Real>::ClassID(); }
  inline static const char* RealName() { return BaseMatlabTraits<_Real>::RealName(); }
};

// Version for char
template<>
struct BaseMatlabTraits<char> : GenericBaseMatlabTraits<char>
{
  enum {
    IsNumeric = 0
  };
  inline static mxClassID ClassID() { return mxCHAR_CLASS; }
  inline static const char* RealName() { return "char"; }
};

// Version for bool
template<>
struct BaseMatlabTraits<bool> : GenericBaseMatlabTraits<bool>
{
  enum {
    IsNumeric = 1
  };
  inline static mxClassID ClassID() { return mxLOGICAL_CLASS; }
  inline static const char* RealName() { return "logical"; }
};


#endif  // DOXYGEN_SHOULD_SKIP_THIS



template <typename T>
struct MatlabTraits : BaseMatlabTraits<T>
{
  inline static mxComplexity Complexity() {
    if (BaseMatlabTraits<T>::IsComplex)
      return mxCOMPLEX;
    else
      return mxREAL;
  }
  inline static const char* ComplexName() {
    if (BaseMatlabTraits<T>::IsNumeric)
      if (BaseMatlabTraits<T>::IsComplex)
	return "complex";
      else
	return "real";
    else  // not numeric
      return "";
  }
  // Type checking
  inline static bool ScalarTypeCheck(const mxArray *pA) {
    return (mxGetClassID(pA) == BaseMatlabTraits<T>::ClassID());
  }
  inline static bool IsA(const mxArray *pA) {
    return ((mxIsComplex(pA) == BaseMatlabTraits<T>::IsComplex) &&
	    ScalarTypeCheck(pA));
  }
  // Allocation
  inline static mxArray* mxAllocator(mwSize ndim,const mwSize *dims) {
    if (BaseMatlabTraits<T>::IsNumeric)
      return mxCreateNumericArray(ndim,dims,BaseMatlabTraits<T>::ClassID(),Complexity());
    else
      return mxCreateCharArray(ndim,dims);
  }
};

#endif  // MatlabTraits_h
