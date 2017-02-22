#!/usr/bin/env python

import argparse
import os
import glob
import subprocess
import errno
import numpy

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

  cmssw_dir = os.path.join(os.environ["CMSSW_BASE"],"src")

  if num_jobs < 1:
    num_jobs = 1

  input_files = [ fullPath(f) for f in glob.glob(os.path.join(input_dir, "*_xsecNom.root")) ]
  num_files = len(input_files)
  input_files = numpy.array_split(numpy.array(input_files), num_jobs)

  num_submitted = 0

  for sublist in input_files:
    if len(sublist) == 0.:
      continue
    job_files = sublist.tolist()
    run_path = os.path.join(run_dir,"scan_point_{}.sh".format(num_submitted))
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
        out_file = os.path.join(output_dir, "limits_and_significances_{}_{}.txt".format(num_submitted, ifile))
        cmd = "./run/scan_point.exe -s -f {} >> {}\n\n".format(f, out_file)
        run_file.write("echo Starting to process file {} of {}\n".format(ifile+1, len(job_files)))
        run_file.write(cmd)

    subprocess.check_call(["JobSubmit.csh",run_path])
    num_submitted += 1

  print("\nSubmitted {} files in {} jobs. Output will be sent to {}.\n".format(
      num_files, num_submitted, output_dir))

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
