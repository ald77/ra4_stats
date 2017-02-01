#!/usr/bin/env python

import argparse
import os
import errno
import glob
import math
import subprocess

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

def SendSignalWorkspaces(input_dir, output_dir, num_jobs, injection_strength, injection_model):
  input_dir = fullPath(input_dir)
  output_dir = fullPath(output_dir)

  run_dir = os.path.join(output_dir, "run")
  ensureDir(run_dir)

  input_files = [ fullPath(f) for f in glob.glob(os.path.join(input_dir, "*SMS*")) ]
  num_files = len(input_files)

  if num_jobs > num_files:
    num_jobs = num_files

  files_per_job = int(math.ceil( float(num_files)/num_jobs ))

  subprocess.check_call(["JobSetup.csh"])

  for ijob in xrange(0,num_jobs):
    run_path = os.path.join(run_dir,"wspace_sig_"+str(ijob)+".sh")
    with open(run_path, "w") as run_file:
      os.fchmod(run_file.fileno(), 0755)
      run_file.write("#! /bin/bash\n\n")

      start_file = ijob*files_per_job
      end_file = min((ijob+1)*files_per_job, num_files)
      for ifile in xrange(start_file, end_file):
        cmd = "./run/wspace_sig.exe -f "+input_files[ifile]+" -o "+output_dir
        cmd += " --sig_strength "+str(injection_strength)+" -u all -l 36.8 -c 0.8484"
        if injection_strength >= 0.:
          cmd += " --unblind none"
          if injection_model != "":
            cmd += " --inject "+injection_model
        run_file.write(cmd+"\n")

    subprocess.check_call(["JobSubmit.csh","./run/wrapper.sh",run_path])

  print "\nSubmitted "+str(num_files)+" files in "+str(num_jobs)+" jobs. Output will be sent to "+output_dir+".\n"

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description = "Submits batch jobs to produce workspaces for each signal mass point",
                                   formatter_class = argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument("output_dir", nargs="?", default = "/net/cms2/cms2r0/babymaker/wspaces/2016_08_10/T1tttt",
                      help="Directory in which to store output workspaces.")
  parser.add_argument("input_dir", nargs="?", default = "/net/cms29/cms29r0/babymaker/babies/2016_08_10/T1tttt/skim_abcd",
                      help="Directory containing input ntuples.")
  parser.add_argument("--num_jobs","-n", type=int, default=50,
                      help="Maximum number of jobs into which to split the mass points")
  parser.add_argument("--injection_strength", type=float, default=-1.,
                      help="Amount of signal to inject. Negative values turn off signal injection. Note that signal injection replaces the data with MC yields, even at injection strength of 0.")
  parser.add_argument("--injection_model", default="",
                      help="Path to signal ntuple to use for signal injection. If unspecified, uses the same signal model as used to construct the likelihood function")
  args = parser.parse_args()

  SendSignalWorkspaces(args.input_dir, args.output_dir, args.num_jobs, args.injection_strength, args.injection_model)
