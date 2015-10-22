#! /bin/bash

if (( "$#" > 0 ))
then
    rootfile=$1
else
    printf '\nSpecify root file with workspace \n\n'
    exit
fi
tempfile=temp.log

combine -M ProfileLikelihood --significance --expectSignal=1 -t -1 $rootfile > $tempfile
echo
echo $rootfile
grep --color=always Significance $tempfile

combine -M Asymptotic -t -1  $rootfile > $tempfile
grep --color=always "Observed Limit" $tempfile
echo

rm $tempfile
