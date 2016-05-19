#!/usr/bin/env python

import os, sys, subprocess
import pprint
import glob
import json
import string
import time


base_folder_name = "variations_05_10_redo/"

resultfolder = base_folder_name+"results/"
model = "T1tttt"

infolder  = "/net/cms27/cms27r0/babymaker/2016_04_29/mc/T1tttt/skim_baseline/"
#fullsim: "/net/cms27/cms27r0/babymaker/2016_04_29/mc/merged_abcd/" 
outfolder = base_folder_name+"wspaces/"
runfolder = base_folder_name+"run/"

PATH = os.getcwd()
if not os.path.isdir(base_folder_name):
    os.mkdir(base_folder_name)
if not os.path.isdir(outfolder):
    os.mkdir(outfolder)
if not os.path.isdir(runfolder):
    os.mkdir(runfolder)
if not os.path.isdir(resultfolder):
    os.mkdir(resultfolder)            

#For all options, need to first ensure there is a systematics file that will be found by wspace_sig.exe
lumis = ["5","7","10","15","20"] 
binnings = ["alternate"]
["nominal","alternate"]


vetoes =["on"]
         #,"off"]
mjthresholds =["350","400"]
mjvariables = ["mj","mj14"]
#["on","off"]
masspoints = ["mGluino-1600_mLSP-1000","mGluino-1400_mLSP-1000","mGluino-1800_mLSP-200"]
#["mGluino-1500_mLSP-100","mGluino-1200_mLSP-800"]
sig_strengths = ["0","1"]

for lumi in lumis: 
    for binning in binnings:
        for mjthres in mjthresholds:
            for mjvar in mjvariables:
                for veto in vetoes:
                    for mass in masspoints:
                        for sig_str in sig_strengths:
                            varname = "lumi"+lumi+"ifb__binning_"+binning+"__veto_"+veto+"_MJdef_"+mjvar+"_threshold_"+mjthres+"_"+mass+"_sig_str"+sig_str
                            exename = runfolder+"wspace_sig_"+varname+".sh"
                            fexe = open(exename,"w")
                            os.system("chmod u+x "+exename)
                            fexe.write("#!/bin/bash\n\n")
                            fexe.write("./run/wspace_sig.exe -l "+lumi+" -b "+binning+" -v "+veto+" -s "+mjthres+" -d "+mjvar+" -g "+sig_str +" -f "+infolder+"fullbaby_SMS-T1tttt_"+mass+"_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_RunIISpring15MiniAODv2-FastAsympt25ns_74X_mcRun2_renorm_abcd.root -o "+outfolder+varname+"\n")
                            # fullsim fexe.write("./run/wspace_sig.exe -l "+lumi+" -b "+binning+" -v "+veto+" -g "+sig_str +" -f "+infolder+"mergedbaby__SMS-T1tttt_"+mass+"__abcd_nfiles_1.root -o "+outfolder+varname+"\n")
                            if sig_str == "1":
                                fexe.write("combine -M ProfileLikelihood --significance --expectSignal=1 "+outfolder+varname+"/wspace_"+model+"_"+mass+"_xsecNom.root | grep Significance | awk '{print $2}' | (read string; echo "+varname+ " $string) > "+resultfolder+varname+".txt \n")
                            if sig_str == "0":
                                fexe.write("combine -M Asymptotic "+outfolder+varname+"/wspace_"+model+"_"+mass+"_xsecNom.root | grep Limit | awk '{print $5}' | (read string; echo "+varname+ " $string) > "+resultfolder+varname+".txt \n")    
                            fexe.close()
                            cmd = "JobSubmit.csh ./run/wrapper.sh "+exename
                            #print cmd
                            os.system(cmd)
