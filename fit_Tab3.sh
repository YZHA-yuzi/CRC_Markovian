#!/bin/bash

#-------------------------------------
# qsub .sh file
#-------------------------------------
for((i=1;i<=2;++i))
do
	for ((b=1;b<=1;++b))
	do 
	Rscript SIM_fitting_Tab3.R "$i" "$b"
	done
done
