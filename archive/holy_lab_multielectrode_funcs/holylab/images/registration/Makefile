MATOPSDIR = $(MATLABINSTALLDIR)/bin
CFLAGS = -g -DMAIN
LDFLAGS = -rpath=$(MATOPSDIR)/glnx86
UTILITYDIR = $(MATLABFUNCDIR)/Utilities
MATRIXDIR = $(MATLABFUNCDIR)/matrix
IMDIR = ..

register_logdetpenalty: register_logdetpenalty.cpp $(IMDIR)/imiterators.cxx $(UTILITYDIR)/mat_variables.h $(UTILITYDIR)/mat_variables.c $(IMDIR)/image_utilities.cpp $(IMDIR)/image_utilities.h
	mex -I$(IMDIR) register_logdetpenalty.cpp $(IMDIR)/image_utilities.cpp
	mex -f $(MATOPSDIR)/matopts.sh $(CFLAGS) -I$(UTILITYDIR) register_logdetpenalty.cpp $(IMDIR)/image_utilities.cpp $(UTILITYDIR)/mat_variables.c -o register_logdetpenalty

register_barrierpenalty_composition: register_barrierpenalty_composition.cpp
	mex -I$(IMDIR) -I$(MATRIXDIR) register_barrierpenalty_composition.cpp

register_logdetpenalty_composition: register_logdetpenalty_composition.cpp $(IMDIR)/imiterators.cxx $(UTILITYDIR)/mat_variables.h $(UTILITYDIR)/mat_variables.c $(MATRIXDIR)/qinterp.cpp jacobianu.cpp
	mex -I$(IMDIR) -I$(MATRIXDIR) register_logdetpenalty_composition.cpp 
#	mex -f $(MATOPSDIR)/matopts.sh $(CFLAGS) -I$(UTILITYDIR) -I$(MATRIXDIR) -I$(IMDIR) register_logdetpenalty_composition.cpp $(UTILITYDIR)/mat_variables.c -o register_logdetpenalty_composition

register_gradE: register_gradE.cpp $(IMDIR)/mat_variables.c $(IMDIR)/mat_variables.h $(IMDIR)/image_utilities.cpp $(IMDIR)/image_utilities.h
	mex -f $(MATOPSDIR)/matopts.sh $(CFLAGS) register_gradE.cpp $(IMDIR)/mat_variables.c $(IMDIR)/image_utilities.cpp -o register_gradE
	mex -I.. register_gradE.cpp $(IMDIR)/image_utilities.cpp

register_gradEreg: register_gradEreg.cpp $(IMDIR)/mat_variables.c $(IMDIR)/mat_variables.h $(IMDIR)/image_utilities.cpp $(IMDIR)/image_utilities.h
	mex -f $(MATOPSDIR)/matopts.sh $(CFLAGS) register_gradEreg.cpp $(IMDIR)/mat_variables.c $(IMDIR)/image_utilities.cpp -o register_gradEreg
	mex -I.. register_gradEreg.cpp $(IMDIR)/image_utilities.cpp

register_g2detJ: register_g2detJ.cpp $(IMDIR)/image_utilities.cpp $(IMDIR)/image_utilities.h
	mex -I.. register_g2detJ.cpp $(IMDIR)/image_utilities.cpp
