const char* rwudata_fields[] = {
  "psim",
  "psit",
  "g",
  "covariant",
  "sqrt",
  "gradpsit",
  "sigma",
  "output",
  "save_memory",
};
int n_rwudata_fields = 9;


struct rwudata_type {
  // Inputs
  const float *psim;        // moving image
  const float *psit;        // target image
  const float *g[];         // deformation
  int n_dims;               // dimensionality of registration problem
  const int *psim_sz;       // size of psim
  const int *psit_sz;       // size of psit
  const int *g_sz;          // size of one element of g{:} (all identical)
  bool covariant;           // boolean covariance flag
  bool sqrt;                // boolean sqrt flag
  bool save_memory;         // sacrifice a bit of accuracy to save memory
  const float *gradpsit;    // spatial gradient of psit
  const double *sigma;      // smoothing length scales

  // Outputs (& their descriptors)
  char *outputID;           // vector of characters describing output arguments
  int n_argsout;            // number of output arguments
  void *outargs[];          // pointers to outputs (of length n_argsout)
};

void struct2rwudata(const mxArray *matstruct,rwudata_type *rwud,int nrhs)
{
  int i,j,nfields;
  const char *curfieldname;
  int curfieldindx;
  int n_dims_m,n_dims_t,n_dims_sigma,n_dims_g;
  const mxArray *gmtrx;
  const int *g_sz_tmp;

  if (!mxIsStruct(matstruct))
    mexErrMsgTxt("Error: input must be a MATLAB structure");
  if (mxGetNumberOfElements(matstruct) != 1)
    mexErrMsgTxt("Error: input must be have only 1 structure element");

  memset(rwud,0,sizeof(rwudata_type));
  // Supply some defaults
  rwud->covariant = true;
  rwud->sqrt = false;
  rwud->save_memory = true;
  
  nfields = mxGetNumberOfFields(matstruct);
  for (i = 0; i < nfields; i++) {
    curfieldname = mxGetFieldNameByNumber(matstruct,i);
    curfieldindx = strmatch(curfieldname, rwudata_fields, n_rwudata_fields);
    curarg = mxGetFieldByNumber(matstruct,0,i);
    switch (curfieldindx) {
    case -1:
      mexPrintf("Fieldname %s not recognized\n",curfieldname);
      mexErrMsgTxt("matstruct_to_rwudata");
      break;
    case 0:   // psim
      if (!mxIsSingle(curarg) || mxIsComplex(curarg))
	mexErrMsgTxt("psim must be of type real single");
      rwud->psim = (float *)mxGetData(curarg);
      n_dims_m = mxGetNumberOfDimensions(curarg);
      rwud->psim_sz = mxGetDimensions(curarg);
      break;
    case 1:   // psit
      if (!mxIsSingle(curarg) || mxIsComplex(curarg))
	mexErrMsgTxt("psit must be of type real single");
      rwud->psit = (float *)mxGetData(curarg);
      n_dims_t = mxGetNumberOfDimensions(curarg);
      rwud->psit_sz = mxGetDimensions(curarg);
      break;
    case 2:   // g
      if (mxIsCell(curarg)) {
	// Cell array input for g
	n_dims_g = mxGetNumberOfElements(curarg);
	rwud->g = (float *) mxMalloc(n_dims_g*sizeof(float *));
	for (j = 0; j < n_dims_g; j++) {
	  // Check the sizes of g elements
	  gmtrx = mxGetCell(curarg,j);
	  if (!mxIsSingle(gmtrx))
	    mexErrMsgTxt("elements of g must be of type single");
	  if (n_dims_g != mxGetNumberOfDimensions(gmtrx))
	    mexErrMsgTxt("dimensionality of g{} is inconsistent");
	  if (j == 0)
	    rwud->g_sz = mxGetDimensions(gmtrx[0]);
	  else {
	    g_sz_tmp = mxGetDimensions(gmtrx);
	    for (k = 0; k < n_dims_g; k++)
	      if (rwud->g_sz[k] != g_sz_tmp[k])
		mexErrMsgTxt("sizes of g{:} are not all the same");
	  }
	  // Get pointers to data
	  rwud->g[j] = (float *) mxGetData(gmtrx);
	}
      } else
	mexErrMsgTxt("For now, g must be a cell array");
      break;
    case 3:   // covariant
      if (!IsScalar(curarg))
	mexErrMsgTxt("covariant must be a scalar");
      rwud->covariant = (bool) mxGetScalar(curarg);
      break;
    case 4:   // sqrt
      if (!IsScalar(curarg))
	mexErrMsgTxt("sqrt must be a scalar");
      rwud->sqrt = (bool) mxGetScalar(curarg);
      break;
    case 5:   // gradpsit
      if (!mxIsSingle(curarg) || mxIsComplex(curarg))
	mexErrMsgTxt("gradpsit must be of type real single");
      rwud->gradpsit = (float *)mxGetData(curarg);
      break;
    case 6:   // sigma
      if (mxIsComplex(curarg))
	mexErrMsgTxt("sigma must be real");
      n_dims_sigma = mxGetNumberOfElements(curarg);  // ndims of _image_
      rwud->sigma = mxGetPr(curarg);
      break;
    case 7:   // output
      if (!mxIsChar(curarg))
	mexErrMsgTxt("output must be a string array");
      rwud->n_argsout = mxGetNumberOfElements(curarg);
      // There's no need for any outputs that are not collected (other
      // than for the "ans" default output), so don't worry about them
      if (nlhs < 1)
	nlhs = 1;
      if (rwud->n_argsout > nlhs)
	rwud->n_argsout = nlhs;
      rwud->outputID = (char *) mxMalloc((rwud->n_argsout+1) * sizeof(char));
      mxGetString(curarg,rwud->outputID,rwud->n_argsout+1);  // don't check output, truncation is OK
      break;
    case 8:   // save_memory
      if (!IsScalar(curarg))
	mexErrMsgTxt("covariant must be a scalar");
      rwud->save_memory = (bool) mxGetScalar(curarg);
      break;
    default:
      mexErrMsgTxt("struct parsing error: this shouldn't happen");
    }
  }

  // Now we have to do some sanity checking.
  // Do the various dimensionalities agree?
  if (rwud->psim)
    rwud->n_dims = n_dims_m;
  if (rwud->psit)
    if (rwud->n_dims) {
      if (rwud->n_dims != n_dims_t)
	mexErrMsgTxt("The dimensionality of psim does not agree with that for psit");
    } else
      rwud->n_dims = n_dims_t;
  if (rwud->sigma)
    if (rwud->n_dims) {
      if (rwud->n_dims != n_dims_sigma)
	mexErrMsgTxt("The dimensionality does not agree with the length of sigma");
    } else
      rwud->n_dims = n_dims_sigma;
  if (rwud->g)
    if (rwud->n_dims) {
      if (rwud->n_dims != n_dims_g)
	mexErrMsgTxt("The dimensionality of g is inconsistent with other parameters");
    } else
      rwud->n_dims = n_dims_g;
  if (rwud->n_dims == 0)
    mexErrMsgTxt("No image data were supplied!");
  // Does the size of g agree with psit?
  if (rwud->psit_sz && rwud->g_sz)
    for (i = 0; i < rwud->n_dims; i++)
      if (rwud->psit_sz[i] != rwud->g_sz[i])
	mexErrMsgTxt("The size of g is not consistent with the size of psit");
    
}

void rwu_work(rwudata_type *rwud)
{
  int warp_output_index,err_output_index,grad_output_index;

  // Somehow, set the output indexes!
  
  calc_jacobian = ??;
  // Everything needs the warped image (even if it is not stored)
  rwu_iminterp(rwud);
  if (rwud->graderr) {
    // We need to calculate the gradient of psig
  if (calc_jacobian) {
    for (pixelIterator = 0; pixelIterator < image.numpix(); pixelIterator++) {
      // This really has to look like the column iterator
      for dimIndex
	for dimIndex
	  Jtmp[] = ;
      detJtmp[pixelIterator] = fabs();  // Store this, because it's small
      /*
      if (rwud->covariant) {
	if (rwud->sqrt)
	  rwud->psig[pixelIterator] *= sqrt(detJ);
	else
	  rwud->psig[pixelIterator] *= detJ;
      }
      */
      if (rwud->graderr) {
	// OK, the pixel iterator has to progress along the dimension
	// for which the gradient in psig is being taken (in cases where
	// graderr is being calculated), so that one can efficiently use
	// it. Copy it backwards (or perhaps just have a circular buffer
	// of pointers?)
	if (pixelIterator > 0) {
	  // Can't do anything if it's the first in the column, but otherwise we can
	}
	if (pixelIterator == last) {
	  // We have to tidy up the last in the column
	}
      }
  } else {
    

  }
