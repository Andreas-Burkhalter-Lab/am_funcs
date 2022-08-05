MATLABFUNCDIR = /home/tim/matlab/matlabfunc
UTILITYDIR = $(MATLABFUNCDIR)/Utilities
MATOPSDIR = /usr/local/matlab2008a/bin
CFLAGS = -g -O -DMAIN -I$(UTILITYDIR) -I/usr/include/valgrind

mzxml_decode: mzxml_decode.cpp $(UTILITYDIR)/mat_variables.h $(UTILITYDIR)/mat_variables.c
	mex mzxml_decode.cpp
	mex -f $(MATOPSDIR)/matopts.sh $(CFLAGS) mzxml_decode.cpp $(UTILITYDIR)/mat_variables.c -o mzxml_decode
