#!/bin/bash

IN_DIR=/bdd/DYAMOND/SAM-4km/OUT_2D
TEMP_DIR=/data/bfildier/DYAMOND_REGIONS/tropics/SAM/tmp
OUT_DIR=/data/bfildier/DYAMOND_REGIONS/tropics/SAM/mean

varid=$1 # a 2D variable


[[ ! -d ${TEMP_DIR} ]] && mkdir ${TEMP_DIR}
cd ${TEMP_DIR}
rm -rf

for file in `ls ${IN_DIR}/*${varid}*`; do

	echo "-- ${file##*/}"

	# temporary file
	end=${file##*_}
	timestamp=${end%%.*}
	outfile=tmp_$timestamp.nc

	# average
	ncwa -a lon,lat -d lat,-30.,30. $file ${TEMP_DIR}/${outfile}

done

# define output file
file_ex=`ls ${IN_DIR}/*${varid}* | head -1`
root=${file_ex##*/}
root=${root%_*}
outfile=${OUT_DIR}/${root}_${varid}.2D.nc

# concatenate
ncrcat *.nc ${outfile}

# exit and remove temp folder
cd ..
rm -rf ${TEMP_DIR}
