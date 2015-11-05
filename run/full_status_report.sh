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

#./run/make_workspace.exe --method m135 -u all --nokappa --no_syst &
#./run/make_workspace.exe --method m135 -u all --no_syst &
#./run/make_workspace.exe --method m135 -u all --use_r4 --no_syst &
./run/make_workspace.exe --method m2l  -u all --lumi 1.264 --no_syst &
./run/make_workspace.exe --method m2l  --lumi 3 --no_syst &
./run/make_workspace.exe --method m1bk --lumi 3 --use_r4 &
./run/make_workspace.exe --method m1bk --lumi 3 --use_r4 --sig_strength 1 &
./run/make_workspace.exe --method m1bk_nodilep --lumi 3 --use_r4 &
./run/make_workspace.exe --method m1bk --lumi 3 --use_r4 --no_syst &
./run/make_workspace.exe --method m1bk -u sideband --lumi 1.264 &

wait

for file in $(ls -A *.root)
do
    ./run/extract_yields.exe $file
done

texify .
