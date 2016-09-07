#!/usr/bin/env python

###### Script to send 58 jobs to the batch finding the limit for each signal mass point
import os, sys, subprocess
import pprint
import glob
import json
import string
import time

# Setting folders
model = "T1tttt"
ntu_date = "2016_08_10"

infolder  = "/net/cms2/cms2r0/babymaker/wspaces/"+ntu_date+"/"+model+"/" 
runfolder = "batch_"+model+"/" 
if not os.path.exists(runfolder):
  os.system("mkdir -p "+runfolder)

cwd = os.getcwd()

#input datasets
inputfiles = [i for i in os.listdir(infolder) if "xsecNom" in i]

os.system("JobSetup.csh")
njobs = 50
files_job = (len(inputfiles)+njobs-1)/njobs
ifile = 0
ijob = 0
for file in inputfiles:
  ifile += 1
  # Creating executable
  if ifile % files_job == 1 or files_job == 1:
    ijob += 1
    exename = runfolder+"/find_limit_sig_"+str(ijob)+".sh"
    fexe = open(exename,"w")
    os.system("chmod u+x "+exename)
    fexe.write("#!/bin/bash\n\n")
    fexe.write(". /cvmfs/cms.cern.ch/cmsset_default.sh \n")
    fexe.write("cd ~/cmssw/CMSSW_7_4_14/src/ \n")
    fexe.write("eval `scramv1 runtime -sh` \n")
    fexe.write("cd "+cwd+" ; \n\n")
  fexe.write("./run/scan_point.exe -f "+infolder+'/'+file+' >> txt/limits_'+model+'_'+str(ijob)+'.txt\n')
  if ifile % files_job == 0 or ifile == len(inputfiles): 
    fexe.close()
    cmd = "JobSubmit.csh ./run/wrapper.sh ./"+exename
    #print cmd
    os.system(cmd)

print "\nSubmitted "+str(ifile)+" files in "+str(ijob)+" jobs\n"
sys.exit(0)
