#!/bin/bash

# Compute low level shear, named LLSU and LLSV for x and y components and LLS for the norm,
# as (varid(3km)-varid(1km))/dz, with approximate heights found in z coordinate

overwrite=False

# Set the directory where your .nc files are located
input_dir="/bdd/DYAMOND/SAM-4km/OUT_3D"
output_dir="/data/bfildier/DYAMOND_REGIONS/tropics/SAM/diagnostics_2D"

i1=16 # = 1078m in SAM's vertical coordinate
i2=27 # = 3175m in SAM's vertical coordinate

# Create a temporary directory to store intermediate output files
tmp_dir=$(mktemp -d)
echo temp dir: $tmp_dir

# show preamble
nchar=48
col1='file processed'
col2='---- elapsed time'
printf "%$(expr ${nchar} + ${#col2})s\n" | tr ' ' '-'
nhalf=$(expr ${nchar} / 2 - ${#col1} / 2)
printf "%${nhalf}s $col1 %${nhalf}s$col2\n"
printf "%$(expr ${nchar} + ${#col2})s\n" | tr ' ' '-'

# missing files
missing_stamps=`comm -23 <(ls ${input_dir}/*_U.nc | awk -F/ '{print $NF}' | cut -f6 -d'_') <(ls ${input_dir}/*_V.nc | awk -F/ '{print $NF}' | cut -f6 -d'_')`
echo missing time stamps: $missing_stamps

# init to keep track of time
SECONDS=0

function secondsToTimeString() 
{
    echo "$(($SECONDS/3600))h:$(($SECONDS%3600/60))m:$(($SECONDS%60))s"
}

# Loop through each netCDF file in the input directory
for file in ${input_dir}/*_U.nc
do

    filename=$(basename $file .nc)
    root=$(basename $file U.nc)
    echo " $filename --- $(secondsToTimeString)" 
    
    # skip if corresponding V file is missing
    time_stamp=`echo $filename | cut -f6 -d'_'`
    [[ ${missing_stamps} =~ .*${time_stamp}.* ]] && echo 'missing V file ; continue' && continue
    
    # Construct the output file path by appending the file name to the temporary directory
    filenameU=${input_dir}/${root}U.nc
    filenameV=${input_dir}/${root}V.nc
    outfilenameU1="${tmp_dir}/${root}U1.nc"
    outfilenameU2="${tmp_dir}/${root}U2.nc"
    outfilenameV1="${tmp_dir}/${root}V1.nc"
    outfilenameV2="${tmp_dir}/${root}V2.nc"
    outfilenamenorm1="${tmp_dir}/${root}norm1.nc"
    outfilenamenorm2="${tmp_dir}/${root}norm2.nc"
    outfilenameUdiff="${tmp_dir}/${root}Udiff.nc"
    outfilenameVdiff="${tmp_dir}/${root}Vdiff.nc"
    outfilenamenormdiff="${tmp_dir}/${root}normdiff.nc"
    
    temp_file="${tmp_dir}/temp.nc"
    output_file_U="${tmp_dir}/${root}LLSU.nc"
    output_file_V="${tmp_dir}/${root}LLSV.nc"
    output_file_norm="${tmp_dir}/${root}LLS.nc"
    
    # skip if output is present
    [[ -f ${output_dir}/$(basename ${output_file_U}) ]] && [[ ! $overwrite == True ]] && echo 'already computed ; continue' && continue
    
    # Use ncks to extract the 2D plan and save it to the temporary directory
    ncks -O -d z,${i1} ${filenameU} ${outfilenameU1}
    ncks -O -d z,${i2} ${filenameU} ${outfilenameU2}
    ncks -O -d z,${i1} ${filenameV} ${outfilenameV1}
    ncks -O -d z,${i2} ${filenameV} ${outfilenameV2}
    
    # compute norm at each level
    # level 1
    ncks -A -v V ${outfilenameV1} ${outfilenameU1}
    ncap2 -s "Unorm=sqrt(U*U+V*V)" ${outfilenameU1} $outfilenamenorm1
    # level 2
    ncks -A -v V ${outfilenameV2} ${outfilenameU2}
    ncap2 -s "Unorm=sqrt(U*U+V*V)" ${outfilenameU2} $outfilenamenorm2
    
    z1=$(ncks -v z ${outfilenameU1} | grep z | tail -1 | grep -oE '[0-9]+\.[0-9]+')
    z2=$(ncks -v z ${outfilenameU2} | grep z | tail -1 | grep -oE '[0-9]+\.[0-9]+')    
    dz=$(echo "$z2 - $z1" | bc)
    
    # difference to get shear
    ncdiff -v U ${outfilenameU1} ${outfilenameU2} ${outfilenameUdiff}
    ncdiff -v V ${outfilenameV1} ${outfilenameV2} ${outfilenameVdiff}
    ncdiff -v Unorm ${outfilenamenorm1} ${outfilenamenorm2} ${outfilenamenormdiff}
    
    ncap2 -s U=U/$dz $outfilenameUdiff $output_file_U
    ncap2 -s V=V/$dz $outfilenameVdiff $output_file_V
    ncap2 -s Unorm=Unorm/$dz $outfilenamenormdiff $output_file_norm
    
    # rename variable
    ncrename -v U,LLSU $output_file_U $temp_file
    mv $temp_file $output_file_U
    ncrename -v V,LLSV $output_file_V $temp_file
    mv $temp_file $output_file_V
    ncrename -v Unorm,LLS $output_file_norm $temp_file
    mv $temp_file $output_file_norm
    
    # rename attributes
    ncatted -O -a long_name,LLSU,o,c,"Low-level shear in U direction (3km-1km)" $output_file_U
    ncatted -O -a long_name,LLSV,o,c,"Low-level shear in V direction (3km-1km)" $output_file_V
    ncatted -O -a long_name,LLS,o,c,"Low-level shear (wind speed at 3km-1km)" $output_file_norm
    ncatted -O -a units,LLSU,o,c,"s-1" $output_file_U
    ncatted -O -a units,LLSV,o,c,"s-1" $output_file_V
    ncatted -O -a units,LLS,o,c,"s-1" $output_file_norm
    ncatted -O -h -a history,global,o,c,"$(date): Shear computed from ${file} between z index 16 (1078m) and z index 27 (3175m) with script computeLowLevelShear ;" $output_file_U
    ncatted -O -h -a history,global,o,c,"$(date): Shear computed from ${file} between z index 16 (1078m) and z index 27 (3175m) with script computeLowLevelShear ;" $output_file_V
    ncatted -O -h -a history,global,o,c,"$(date): Shear computed from ${file} between z index 16 (1078m) and z index 27 (3175m) with script computeLowLevelShear ;" $output_file_norm
    
    # remove temporary files
    rm $outfilenameU1 $outfilenameU2 $outfilenameV1 $outfilenameV2 $outfilenamenorm1 $outfilenamenorm2 $outfilenameUdiff $outfilenameVdiff $outfilenamenormdiff
    
    mv $tmp_dir/* $output_dir/
    
done


# Remove the temporary directory and its contents
rm -r "${tmp_dir}"
