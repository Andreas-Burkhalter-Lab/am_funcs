UTILITYDIR = $(MATLABFUNCDIR)/Utilities
IMAGEDIR = $(MATLABFUNCDIR)/images
MATOPSDIR = $(MATLABINSTALLDIR)/bin
CFLAGS = -g -DMAIN -O

bn_image_filter: bn_image_filter.cpp $(IMAGEDIR)/imiterators.cxx $(UTILITYDIR)/mat_variables.h $(UTILITYDIR)/mat_variables.c #$(UTILITYDIR)/newdelete.h $(UTILITYDIR)/newdelete.cpp
	mex -msse2 -I$(UTILITYDIR) -I$(IMAGEDIR) bn_image_filter.cpp $(UTILITYDIR)/newdelete.cpp
	mex -f $(MATOPSDIR)/matopts.sh -msse2 $(CFLAGS) -I$(UTILITYDIR) -I$(IMAGEDIR) bn_image_filter.cpp $(UTILITYDIR)/mat_variables.c -o bn_image_filter

msams_converge1: msams_converge1.cpp msams_core.cpp $(UTILITYDIR)/mat_variables.h $(UTILITYDIR)/mat_variables.c $(UTILITYDIR)/newdelete.h $(UTILITYDIR)/newdelete.cpp
	mex -I$(UTILITYDIR) msams_converge1.cpp $(UTILITYDIR)/newdelete.cpp
	mex -f $(MATOPSDIR)/matopts.sh $(CFLAGS) -I$(UTILITYDIR) msams_converge1.cpp $(UTILITYDIR)/newdelete.cpp $(UTILITYDIR)/mat_variables.c -o msams_converge1

msams_converge_mex: msams_converge_core.h msams_converge_core.cpp msams_converge_mex.cpp landmarked_neighbors.cpp landmarked_neighbors.h $(UTILITYDIR)/mat_variables.h $(UTILITYDIR)/mat_variables.c
	mex msams_converge_mex.cpp
	mex -f $(MATOPSDIR)/matopts.sh $(CFLAGS) msams_converge_mex.cpp $(UTILITYDIR)/mat_variables.c -o msams_converge
#	mex $(CFLAGS) msams_converge_mex.cpp $(UTILITYDIR)/mat_variables.c -o msams_converge

#all: msams_converge_mex
