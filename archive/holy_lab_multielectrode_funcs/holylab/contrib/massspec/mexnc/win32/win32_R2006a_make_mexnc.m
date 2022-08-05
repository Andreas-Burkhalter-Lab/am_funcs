cd ..
mex -v -f win/win32_R2006a.bat -output mexnc mexgateway.c netcdf2.c netcdf3.c common.c
cd win32
