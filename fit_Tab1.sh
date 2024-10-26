#!/bin/bash

#-------------------------------------
# qsub .sh file
#-------------------------------------
for((s=1;s<=2;++s))
do
	for((i=1;i<=2;++i))
	do
		for ((b=1;b<=2;++b))
		do 
		Rscript SIM_fitting_Tab1.R "$s" "$i" "$b"
		done
	done
done

