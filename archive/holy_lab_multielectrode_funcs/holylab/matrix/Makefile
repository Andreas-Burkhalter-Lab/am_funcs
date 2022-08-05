UTILITYDIR = $(MATLABFUNCDIR)/Utilities
IMAGEDIR = $(MATLABFUNCDIR)/images
MATOPSDIR = $(MATLABINSTALLDIR)/bin
CFLAGS = -g -DMAIN -O

rotate_quaternion: rotate_quaternion.cpp $(UTILITYDIR)/mat_variables.h $(UTILITYDIR)/mat_variables.c 
	mex -O rotate_quaternion.cpp
	mex -f $(MATOPSDIR)/matopts.sh $(CFLAGS) -I$(UTILITYDIR) rotate_quaternion.cpp $(UTILITYDIR)/mat_variables.c -o rotate_quaternion

cummax: cummax.cpp $(IMAGEDIR)/imiterators.cxx $(UTILITYDIR)/mat_variables.h $(UTILITYDIR)/mat_variables.c
	mex -O -I$(IMAGEDIR) cummax.cpp
