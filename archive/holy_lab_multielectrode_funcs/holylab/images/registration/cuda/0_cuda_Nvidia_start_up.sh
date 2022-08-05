#!/bin/bash 

# bash script to initiaze Nvidia GPU device if not initized by default
# "ls /dev/nvidia*" will display the number of GPU cards installed
# if not all GPU devices show, execute this bash script using sudo bash

/sbin/modprobe nvidia 
if [ "$?" -eq 0 ]; then 
	# Count the number of NVIDIA controllers found. 
	NVDEVS=`lspci | grep -i NVIDIA` 
	N3D=`echo "$NVDEVS" | grep "3D controller" | wc -l` 
	NVGA=`echo "$NVDEVS" | grep "VGA compatible controller" | wc -l` 
	N=`expr $N3D + $NVGA - 1` 
	for i in `seq 0 $N`; do 
		mknod -m 666 /dev/nvidia$i c 195 $i 
	done 

	mknod -m 666 /dev/nvidiactl c 195 255 
else 
	exit 1 
fi
