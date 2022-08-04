# Tims standard Windows NMAKE file, extensively modified to provide
# for a build and update of the SON filing system as required for Windows
#
# Usage:     NMAKE option -f son16.mak
#
# option:    DEBUG=[0|1]  (DEBUG not defined is equivalent to DEBUG=1)
#
# First of all, define variables to control C compilation and linking

!if "$(DEBUG)" != "0"
CFLAGS = -DDEBUG -DSTRICT -Od -Zip -AMw -G2sw -W4 -Gc -c
LFLAGS = /NOE/NOD/CODEVIEW
!ELSE
CFLAGS = -AMw -W4 -G2sw -Oselgit -Zp -c
LFLAGS = /NOE/NOD
!ENDIF

# The primary target is son.dll, so we put this one first

UPDATE: son16.dll

son16.dll : sondll.obj son16.def son16.mak
    link $(LFLAGS) sondll libentry, son16.dll,,mdllcew libw, son16.def
    rc son16.dll
    implib son16.lib son16.def

# Rebuild the object files as necessary

sondll.obj : sondll.c son.c son.h sonintl.h son16.mak
    cl $(CFLAGS) sondll.c >son.err
    type son.err

