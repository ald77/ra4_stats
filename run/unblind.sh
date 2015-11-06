#! /bin/bash

shopt -s nullglob

texify(){
    orig_dir=$(pwd)
    cd $1
    for i in $(ls -A)
    do
	if [ -f $i ] && [ "${i: -4}" == ".tex" ]
	then
            pdfname=$(basename "$i" .tex).pdf
	    if [ ! -f $pdfname ] || [ $i -nt $pdfname ]
	    then
		echo "Compiling $i"
		pdflatex --shell-escape $i 1> /dev/null
		auxname=$(basename "$i" .tex).aux
		logname=$(basename "$i" .tex).log
		rm -f $auxname $logname
	    else
		echo "Skipping $i"
	    fi
	fi
    done
    cd $orig_dir
}

./compile.sh

./run/make_workspace.exe --method m1bk -l 1.264 -u sideband &
./run/make_workspace.exe --method m1bk -l 1.264 -u 1b &
./run/make_workspace.exe --method m1bk -l 1.264 -u all &
./run/make_workspace.exe --method m1bk -l 1.264 -u all --use_r4 &

wait

for file in $(ls -A *.root)
do
    ./run/extract_yields.exe -f $file
done

texify .
