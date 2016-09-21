#!/usr/bin/env python

import argparse
import os
import glob
import subprocess
import errno
import math

def ensureDir(path):
  try:
    os.makedirs(path)
  except OSError as e:
    if e.errno == errno.EEXIST and os.path.isdir(path):
      pass
    else:
      raise

def fullPath(path):
  return os.path.realpath(os.path.abspath(os.path.expanduser(path)))

def SendLimits(input_dir, output_dir, num_jobs):
  input_dir = fullPath(input_dir)
  output_dir = fullPath(output_dir)

  run_dir = os.path.join(output_dir, "run")
  ensureDir(run_dir)

  input_files = [ fullPath(f) for f in glob.glob(os.path.join(input_dir, "*_xsecNom.root")) ]
  num_files = len(input_files)

  if num_jobs > num_files:
    num_jobs = num_files

  files_per_job = int(math.ceil( float(num_files)/num_jobs ))

  subprocess.check_call(["JobSetup.csh"])
  cwd = os.getcwd()

  for ijob in xrange(0,num_jobs):
    run_path = os.path.join(run_dir,"find_limit_sig_"+str(ijob)+".sh")
    with open(run_path, "w") as run_file:
      os.fchmod(run_file.fileno(), 0755)
      run_file.write("#! /bin/bash\n\n")
      run_file.write(". /cvmfs/cms.cern.ch/cmsset_default.sh\n")
      run_file.write("cd ~/cmssw/CMSSW_7_4_14/src\n")
      run_file.write("eval `scramv1 runtime -sh`\n")
      run_file.write("cd "+cwd+"\n\n")
      start_file = ijob*files_per_job
      end_file = min((ijob+1)*files_per_job, num_files)
      for ifile in xrange(start_file, end_file):
        out_file = os.path.join(output_dir, "limits_and_significances_"+str(ifile)+".txt")
        cmd = "./run/scan_point.exe -s -f "+input_files[ifile]+" >> "+out_file
        run_file.write(cmd+"\n")

    subprocess.check_call(["JobSubmit.csh","./run/wrapper.sh",run_path])

  print "\nSubmitted "+str(num_files)+" files in "+str(num_jobs)+" jobs. Output will be sent to "+output_dir+".\n"


if __name__ == "__main__":
  parser = argparse.ArgumentParser(description="Submits batch jobs to compute limits and significances from existing workspaces",
                                   formatter_class = argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument("output_dir", nargs="?", default = "txt",
                      help = "Directory in which to store computed results.")
  parser.add_argument("input_dir", nargs="?", default = "/net/cms2/cms2r0/babymaker/wspaces/2016_08_10/T1tttt",
                      help = "Directory containing workspaces to be processed.")
  parser.add_argument("--num_jobs","-n", type=int, default=50,
                      help = "nNumber of jobs into which to split processing of workspaces")
  args = parser.parse_args()

  SendLimits(args.input_dir, args.output_dir, args.num_jobs)
