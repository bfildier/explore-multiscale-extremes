#!/bin/bash

for file in `ls logs/run*.o`
do
	echo 
	echo $file
       	echo
       	head -19 $file
	echo '*'
	echo '*'
	tail -3 $file
       	echo
       	echo ------------------
done

