UTILITYDIR = $(MATLABFUNCDIR)/Utilities
IMUTILITYDIR = $(MATLABFUNCDIR)/images
MATOPSDIR = $(MATLABINSTALLDIR)/bin
CFLAGS = -g -DMAIN -I$(UTILITYDIR) -I$(IMUTILITYDIR) -I/usr/include/valgrind

array_restrict: array_restrict.cpp $(IMUTILITYDIR)/imiterators.cxx $(IMUTILITYDIR)/image_utilities.cpp $(UTILITYDIR)/mat_variables.h $(UTILITYDIR)/mat_variables.c
	mex -I$(IMUTILITYDIR) array_restrict.cpp $(IMUTILITYDIR)/image_utilities.cpp
	mex -f $(MATOPSDIR)/matopts.sh $(CFLAGS) array_restrict.cpp $(IMUTILITYDIR)/image_utilities.cpp $(UTILITYDIR)/mat_variables.c -o array_restrict


array_prolong: array_prolong.cpp $(IMUTILITYDIR)/imiterators.cxx $(IMUTILITYDIR)/image_utilities.cpp $(UTILITYDIR)/mat_variables.h $(UTILITYDIR)/mat_variables.c $(UTILITYDIR)/newdelete.cpp
	mex -I$(IMUTILITYDIR) -I$(UTILITYDIR) array_prolong.cpp $(IMUTILITYDIR)/image_utilities.cpp $(UTILITYDIR)/newdelete.cpp
	mex -f $(MATOPSDIR)/matopts.sh $(CFLAGS) array_prolong.cpp $(IMUTILITYDIR)/image_utilities.cpp $(UTILITYDIR)/mat_variables.c -o array_prolong


array_quadratic_penalty: array_quadratic_penalty.cpp $(IMUTILITYDIR)/imiterators.cxx $(UTILITYDIR)/MatlabIO.h $(UTILITYDIR)/MatlabTraits.h
	mex -I$(IMUTILITYDIR) -I$(UTILITYDIR) array_quadratic_penalty.cpp
#	mex -f $(MATOPSDIR)/matopts.sh $(CFLAGS) array_quadratic_penalty.cpp -o array_quadratic_penalty
