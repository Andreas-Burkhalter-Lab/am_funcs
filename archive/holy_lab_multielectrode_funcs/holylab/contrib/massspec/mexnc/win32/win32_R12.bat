@echo off
rem MSVC70_R12.BAT
rem
rem    Compile and link options used for building mexnc.dll 
rem    using the Microsoft Visual C++ compiler version 7.1 and Matlab 6.0
rem    and single threaded NetCDF 3.6.2
rem    
rem
********************************************************************
set UNIDATA_INC=h:\src\netcdf-3.6.2-snapshot2007072001\libsrc
set UNIDATA_LIB=h:\src\netcdf-3.6.2-snapshot2007072001\win32\NET\debug
set UNIDATA_LIBS=h:\src\netcdf-3.6.2-snapshot2007072001\win32\NET\debug\netcdf.lib
rem ********************************************************************
rem General parameters
rem ********************************************************************

set MATLAB=%MATLAB%
set MSVCDir=%MSVCDir%
set MSDevDir=%MSVCDir%
set PATH=%MSDevDir%\bin;%PATH%
set INCLUDE=%UNIDATA_INC%;%MSVCDir%\INCLUDE;%INCLUDE%
set LIB=%MSVCDir%\LIB;%MSVCDir%\PlatformSDK\LIB;%LIB%

rem ********************************************************************
rem Compiler parameters
rem ********************************************************************
set COMPILER=cl
set COMPFLAGS= -DMATLAB_MEX_FILE -DDLL_NETCDF -DMEXNCR12 /nologo /W3 /GR /GX /c
set OPTIMFLAGS=/ML 
set DEBUGFLAGS=/MDd -Zi -Fd"%OUTDIR%%MEX_NAME%.pdb"
set NAME_OBJECT=/Fo

rem ********************************************************************
rem Linker parameters
rem ********************************************************************
set LIBLOC=%MATLAB%\extern\lib\win32\microsoft\msvc60
set LINKER=link
set LINKFLAGS=/dll /NODEFAULTLIB:msvcrt.lib /VERBOSE:lib /export:%ENTRYPOINT% /MAP /LIBPATH:"%UNIDATA_LIB%" %UNIDATA_LIBS% /LIBPATH:"%LIBLOC%" libmx.lib libmex.lib libmat.lib /implib:%LIB_NAME%.x
set LINKOPTIMFLAGS=
set LINKDEBUGFLAGS=/debug
set LINK_FILE=
set LINK_LIB=
set NAME_OUTPUT=/out:"%OUTDIR%%MEX_NAME%.dll"
set RSP_FILE_INDICATOR=@

rem ********************************************************************
rem Resource compiler parameters
rem ********************************************************************
set RC_COMPILER=rc /fo "%OUTDIR%mexversion.res"
set RC_LINKER=

set POSTLINK_CMDS=del "%OUTDIR%%MEX_NAME%.map"
set POSTLINK_CMDS1=del %LIB_NAME%.x
