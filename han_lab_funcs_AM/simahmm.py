import sys

# def markov(fname, nplane):
import sima
import sima.motion
import sima.segment
import numpy as np
# import matplotlib.pyplot as plt

# return fname
 #File location here
 #%cd E:\2016
 #Filename here
 #name='2016160118_000_001'
#Number of planes here
# #nplane=4

fname=sys.argv[1]
#seq=[sima.Sequence.create('HDF5',fname+'.mat','txy') for i in range(1,nplane+1)]
print(fname+'.mat')
print([fname[:end-7]+'hmm'+fname[end-4]+'.hdf5'])
# #HMM transformation, more accurate
#mchmm = sima.motion.HiddenMarkov2D(granularity='row',max_displacement=[20, 30], verbose=False)
#dshmm = mchmm.correct(seq, 'hmm.sima')
#dshmm.export_frames([fname[:end-7]+'hmm'+str(i)+'.hdf5' for i in range(1,5)],fmt='HDF5')