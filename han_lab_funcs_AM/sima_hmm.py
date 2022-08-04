import sys

print(str(sys.argv))


# #File location here
# %cd E:\2016
# #Filename here
# name='2016160118_000_001'
# #Number of planes here
# nplane=4

# import sima
# import sima.motion
# import sima.segment
# import numpy as np
# import matplotlib.pyplot as plt
# %matplotlib inline
# seq=[sima.Sequence.create('HDF5',name+'_plane'+str(i)+'.mat','txy') for i in range(1,nplane+1)]

# #HMM transformation, more accurate
# mchmm = sima.motion.HiddenMarkov2D(granularity='row',max_displacement=[20, 30], verbose=False)
# dshmm = mchmm.correct(seq, 'hmm.sima')
# dshmm.export_frames([name+'hmm'+str(i)+'.hdf5' for i in range(1,5)],fmt='HDF5')