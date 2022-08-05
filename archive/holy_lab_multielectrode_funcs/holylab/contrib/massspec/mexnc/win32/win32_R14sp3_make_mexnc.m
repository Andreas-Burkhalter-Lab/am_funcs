cd ..
mex -v -f win/win32_R14sp3.bat -output mexnc mexgateway.c netcdf2.c netcdf3.c common.c
cd win32
