/* The SONDLL source simply provides a very simple wrapper for the DOS      */
/* SON code. This must be compiled with -AM -GD -G2s, using C7 for 16bit    */

#include <windows.h>

#include "son.c"

#ifdef WIN32

/****************************************************************************
   FUNCTION: DllMain(HANDLE, DWORD, LPVOID)

   PURPOSE:  DllMain is called by Windows when
             the DLL is initialized, Thread Attached, and other times.
             Refer to SDK documentation, as to the different ways this
             may be called.
             
             The DllMain function should perform additional initialization
             tasks required by the DLL.  In this example, no initialization
             tasks are required.  DllMain should return a value of 1 if
             the initialization is successful.
           
*******************************************************************************/
INT  APIENTRY DllMain(HANDLE hInst, DWORD ul_reason_being_called, LPVOID lpReserved)
{
   if (ul_reason_being_called == DLL_PROCESS_ATTACH)
      SONInitFiles();

   return 1;

	UNREFERENCED_PARAMETER(hInst);
	UNREFERENCED_PARAMETER(lpReserved);
}

#else

/****************************************************************************
   FUNCTION: LibMain(HANDLE, WORD, WORD, LPSTR)

   PURPOSE:  LibMain is called by Windows when the DLL is loaded and
             initialized.
             
             The LibMain function should perform additional initialization
             tasks required by the DLL.  In this example, initialization
             is calling InitFiles.  LibMain should return a value of 1 if
             the initialization is successful.
           
*******************************************************************************/
int WINAPI LibMain(HANDLE hInst, WORD wDataSeg, WORD wHeapSize,
                        LPSTR lpszCmdLine)
{
/*	UNREFERENCED_PARAMETER(hInst);      What is equivalent for WIN16 ?
	UNREFERENCED_PARAMETER(wDataSeg);
	UNREFERENCED_PARAMETER(wHeapSize);
	UNREFERENCED_PARAMETER(lpszCmdLine); */

    SONInitFiles();

    if (wHeapSize > 0)
        UnlockData(0) ;

    return 1;
}

#endif

void WINAPI WEP( int nParameter)
{
    return;
}


