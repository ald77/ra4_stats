#! /bin/bash

shopt -s nullglob

if [ $# -lt 1 ]
then
    echo "Must specify an input directory"
    exit 1
fi

lim_dir=$1
num_parallels=8
if [ $# -gt 1 ]
then
    num_parallels=$2
fi

index=0
echo "-----" > txt/limits.txt
for file in $(ls -A $lim_dir/*_xsecNom.root)
do
    echo "Processing $file"
    index=$((index+1))
    ./run/scan_point.exe -f $file < /dev/null | tail -n 1 >> txt/limits.txt &
    if (( $index % $num_parallels == 0 )) && [ $index -ne 0 ]
    then
        echo "Waiting for running jobs to finish..."
        wait
    fi
done

echo "Waiting for running jobs to finish..."
wait
