#!/bin/bash

for file in `ls logs/run*.e`
do
	echo 
	echo $file
       	echo
       	tail -15 $file
       	echo
       	echo ------------------
done

