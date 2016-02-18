#! /bin/bash

orig_dir=$(pwd)
mkdir -p ~/cmssw
cd ~/cmssw
scramv1 project CMSSW CMSSW_7_4_14 || mkdir -p ~/cmssw/CMSSW_7_4_14/src
scramv1 project CMSSW CMSSW_7_1_5 || mkdir -p ~/cmssw/CMSSW_7_1_5/src
cd ~/cmssw/CMSSW_7_1_5/src
eval `scramv1 runtime -sh`
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
git checkout v5.0.1
scram b -j $(getconf _NPROCESSORS_ONLN) -k
cd ~/cmssw/CMSSW_7_4_14/src
eval `scramv1 runtime -sh`
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
git checkout 74x-root6
scram b -j $(getconf _NPROCESSORS_ONLN) -k
cd $orig_dir
./compile.sh
