extern int validate_dimensions(const int *,int,const int *,int);
extern int skip_unity_dimensions(const int *din,int n,int *dout,int noutmax);
extern int skip_unity_dimensions_index(const int *din,int n,int *indxout,int noutmax);
extern long calc_pixel_skip(const int *sz,int n,long *pixel_skip);
extern int isVector(const int *sz,const int n_dims);
