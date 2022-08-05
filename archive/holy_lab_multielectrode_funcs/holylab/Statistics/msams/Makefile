UTILITYDIR = $(MATLABFUNCDIR)/Utilities
IMAGEDIR = $(MATLABFUNCDIR)/images
STATSDIR = $(MATLABFUNCDIR)/Statistics
MATOPSDIR = $(MATLABINSTALLDIR)/bin
CFLAGS = -g -DMAIN

msams_converge1: msams_converge1.cpp msams_core.cpp $(UTILITYDIR)/mat_variables.h $(UTILITYDIR)/mat_variables.c $(UTILITYDIR)/newdelete.h $(UTILITYDIR)/newdelete.cpp
	mex -I$(UTILITYDIR) msams_converge1.cpp $(UTILITYDIR)/newdelete.cpp
	mex -f $(MATOPSDIR)/matopts.sh $(CFLAGS) -I$(UTILITYDIR) msams_converge1.cpp $(UTILITYDIR)/newdelete.cpp $(UTILITYDIR)/mat_variables.c -o msams_converge1

msams_converge_mex: msams_converge_core.h msams_converge_core.cpp msams_converge_mex.cpp $(STATSDIR)/landmarked_neighbors.cpp $(STATSDIR)/landmarked_neighbors.h $(UTILITYDIR)/mat_variables.h $(UTILITYDIR)/mat_variables.c
	mex -I$(STATSDIR) msams_converge_mex.cpp
	mex -f $(MATOPSDIR)/matopts.sh $(CFLAGS) -I$(UTILITYDIR) -I$(STATSDIR) msams_converge_mex.cpp $(UTILITYDIR)/mat_variables.c -o msams_converge_mex
#	mex $(CFLAGS) msams_converge_mex.cpp $(UTILITYDIR)/mat_variables.c -o msams_converge

#all: msams_converge_mex
