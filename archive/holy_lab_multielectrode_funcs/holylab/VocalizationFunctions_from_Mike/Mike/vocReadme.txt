This document explains the data structure output in vocalizationRun


data.wavfile			This is the the complete path to the audio file. 
						Run wavread(data.wavefile) to load into memory or for other analysis.
data.sngfile			This is the complete path to the sonogram file.
						Run ReadSonogram(data.sngfile) to load the SNG into memory.
						Run spsngplot(data.sngfile) to plot the sonogram with the slider window.
						Run whistimesplot(data.sngfile,[],data.twhis) to plot call start/stop times with red lines
data.tacq				This is the total acquisition time in seconds for the wav file, from the SNG header.
						Calculate data.nwhis/data.tacq/60 to determine the call rate in calls per minute.

						
data.nwhis				This is the # of USVs in the file, based on whistimes.
data.twhis				This is a 2XN vector of [start;stop] where each column is a USV.
						Duration = data.twhis(2,:) - data.twhis(1,:)
						Pause = data.twhis(1,2:end) - data.twhis(2,1:end-1)

data.nobj				This is the number of spectrogram objects determined for each call.
						
data.objmf				This is the mean frequency in Hz of each spectrogram object, for each call.