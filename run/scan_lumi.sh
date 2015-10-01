#! /bin/bash
./compile.sh 

./run/make_workspace.exe -l 0.5 &
./run/make_workspace.exe -l 1 &
./run/make_workspace.exe -l 1.5 &
./run/make_workspace.exe -l 2 &
./run/make_workspace.exe -l 2.5 &
./run/make_workspace.exe -l 3 &
./run/make_workspace.exe -l 0.5 --no_syst &
./run/make_workspace.exe -l 1 --no_syst &
./run/make_workspace.exe -l 1.5 --no_syst &
./run/make_workspace.exe -l 2 --no_syst &
./run/make_workspace.exe -l 2.5 --no_syst &
./run/make_workspace.exe -l 3 --no_syst &

wait