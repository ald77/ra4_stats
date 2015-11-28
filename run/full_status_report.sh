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

./run/make_workspace.exe --method m2l  -u all --lumi 2.1 --no_syst &
./run/make_workspace.exe --method m1bk --lumi 2.1 --use_r4 &
./run/make_workspace.exe --method m1bk --lumi 2.1 --sig_strength 1 --use_r4 &
./run/make_workspace.exe --method m1bk_nodilep --lumi 2.1 --use_r4 &
./run/make_workspace.exe --method m1bk --lumi 2.1 --no_syst --use_r4 &
./run/make_workspace.exe --method m1bk -u all --lumi 2.1 &
./run/make_workspace.exe --method m1bk -u all --lumi 2.1 --use_r4 &

wait

for file in $(ls -A *.root)
do
#    ./run/extract_yields.exe -f $file
done

texify .
