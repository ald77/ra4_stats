#!/usr/bin/env python

from __future__ import print_function

import argparse
import os
import errno
import glob
import numpy
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

  cmssw_dir = os.path.join(os.environ["CMSSW_BASE"],"src")

  if num_jobs < 1:
    num_jobs = 1

  input_files = [ fullPath(f) for f in glob.glob(os.path.join(input_dir, "*SMS*")) ]
  num_files = len(input_files)
  input_files = numpy.array_split(numpy.array(input_files), num_jobs)

  num_submitted = 0

  for sublist in input_files:
    if len(sublist) == 0:
      continue
    job_files = sublist.tolist()
    run_path = os.path.join(run_dir,"wspace_sig_{}.sh".format(num_submitted))
    with open(run_path, "w") as run_file:
      os.fchmod(run_file.fileno(), 0755)
      run_file.write("#! /bin/bash\n\n")

      run_file.write("DIRECTORY=`pwd`\n")
      run_file.write("cd {}\n".format(cmssw_dir))
      run_file.write(". /net/cms2/cms2r0/babymaker/cmsset_default.sh\n")
      run_file.write("eval `scramv1 runtime -sh`\n")
      run_file.write("cd $DIRECTORY\n\n")

      for ifile in range(len(job_files)):
        f = job_files[ifile]
        cmd = "./run/wspace_sig.exe -f {} -o {} --sig_strength {} -u all -l 35.9 -p".format(
          f, output_dir, (injection_strength if injection_strength >= 0. else 0.))
        if injection_strength >= 0.:
          cmd += " --unblind none"
          if injection_model != "":
            cmd += " --inject "+injection_model
        run_file.write("echo Starting to process file {} of {}\n".format(ifile+1, len(job_files)))
        run_file.write(cmd+"\n\n")

    subprocess.check_call(["JobSubmit.csh",run_path])
    num_submitted += 1

  print("\nSubmitted {} files in {} jobs. Output will be sent to {}.\n".format(
      num_files, num_submitted, output_dir))

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
