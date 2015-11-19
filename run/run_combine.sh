#! /bin/bash

if (( "$#" > 0 ))
then
    rootfile=$1
else
    printf '\nSpecify root file with workspace \n\n'
    exit
fi
tempfile=temp.log

printf '\n'$rootfile'\n\n' 
combine -M ProfileLikelihood --significance --expectSignal=1 $rootfile > $tempfile
echo Observed `grep --color=always Significance $tempfile`

combine -M Asymptotic $rootfile > $tempfile
grep --color=always "Observed Limit" $tempfile
echo A posteriori `grep --color=always "Expected 16.0%" $tempfile`
echo A posteriori `grep --color=always "Expected 50.0%" $tempfile`
echo A posteriori `grep --color=always "Expected 84.0%" $tempfile`
echo

if (( "$#" > 1 ))
then
	combine -M ProfileLikelihood --significance --expectSignal=1  --toysFreq $rootfile > $tempfile
	echo A posteriori expected `grep --color=always Significance $tempfile`
	
	combine -M ProfileLikelihood --significance --expectSignal=1 -t -1  --toysFreq $rootfile > $tempfile
	echo "A posteriori expected (with -t -1)" `grep --color=always Significance $tempfile`
	
	combine -M ProfileLikelihood --significance --expectSignal=1 -t -1 $rootfile > $tempfile
	echo A priori expected `grep --color=always Significance $tempfile`
	
	combine -M Asymptotic -t -1  $rootfile > $tempfile
	echo A priori `grep --color=always "Expected 16.0%" $tempfile`
	echo A priori `grep --color=always "Expected 50.0%" $tempfile`
	echo A priori `grep --color=always "Expected 84.0%" $tempfile`
	echo
fi

rm $tempfile
