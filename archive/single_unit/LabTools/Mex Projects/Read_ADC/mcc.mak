bin_target : main.dll

main.obj : "C:\Program Files\Microsoft Visual Studio\MyProjects\Read_ADC\main.c"
	cl  -c -Zp8 -G5 -W3 -DMATLAB_MEX_FILE -nologo /Fomain.obj -IC:\matlabR12\extern\include  -O2 -Oy- -DNDEBUG "C:\Program Files\Microsoft Visual Studio\MyProjects\Read_ADC\main.c"

main.dll : main.obj
	rc /fo "mexversion.res"  "C:\matlabR12\extern\include\mexversion.rc"
	link /out:"main.dll" /dll /export:mexFunction /MAP /LIBPATH:"C:\matlabR12\extern\lib\win32\microsoft\msvc60" libmx.lib libmex.lib libmatlbmx.lib libmat.lib /implib:_lib7726.x  @main_master.rsp 
	del "main.map"
	del _lib7726.x
	if exist "mexversion.res" del "mexversion.res"

