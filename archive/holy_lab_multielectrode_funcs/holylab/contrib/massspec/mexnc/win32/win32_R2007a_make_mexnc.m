cd ..
mex -v -f win32/win32_R2007a.bat -output mexnc mexgateway.c netcdf2.c netcdf3.c common.c
cd win32
