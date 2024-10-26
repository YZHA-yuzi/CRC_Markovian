#!/bin/bash

#-------------------------------------
# qsub .sh file
#-------------------------------------
for((i=1;i<=3;++i))
do
	Rscript RealApp.R "$i" 
done
