#! /usr/bin/env python

import argparse
import os
import shutil
import re
import tempfile
import subprocess
import ROOT
import sys

def GetRegex(name):
    return re.search("t1tttt_(.*?)_(.*?)_lumi_(.*?)_met_(.*?)_inf_njets_(.*?)_inf_nbm_(.*?)_inf_tkveto_(.*?).root", name)

def GetParams(name):
    param_strings = GetRegex(name)
    mglu = int(param_strings.group(1))
    mlsp = int(param_strings.group(2))
    lumi = float(param_strings.group(3))
    met = float(param_strings.group(4))
    njets = int(param_strings.group(5))
    nbm = int(param_strings.group(6))
    tkveto = param_strings.group(7)
    if tkveto == "true":
        tkveto = True
    else:
        tkveto = False
    return (mglu, mlsp, lumi, met, njets, nbm, tkveto)

def GetLimit(filename):
    f = ROOT.TFile(filename, "read")
    for wp in f.limit:
        if wp.quantileExpected == 0.5:
            return wp.limit

    return -1.

def GetSignificance(filename):
    f = ROOT.TFile(filename, "read")
    for wp in f.limit:
        if wp.quantileExpected == -1.:
            return wp.limit

    return -1.

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extracts binning parameters, limits, and significance from a single (aggregate) bin workspace")
    parser.add_argument("-f", "--file", required=True, help="File from which to extract parameters")
    parser.add_argument("-o", "--output", default="", help="Directory in which to store results")
    args = parser.parse_args()

    params = GetParams(args.file)

    origdir = os.getcwd()
    workdir = tempfile.mkdtemp()
    name = os.path.basename(args.file)
    os.symlink(args.file, os.path.join(workdir, name))
    os.chdir(workdir)
    
    subprocess.call(["combine","-M","Asymptotic","-t","-1","--toysFreq",name])
    subprocess.call(["combine","-M","ProfileLikelihood","--significance","--expectSignal=1","-t","-1","--toysFreq",name])

    limit = GetLimit(os.path.join(workdir, "higgsCombineTest.Asymptotic.mH120.root"))
    signif = GetSignificance(os.path.join(workdir, "higgsCombineTest.ProfileLikelihood.mH120.root"))

    params += (limit, signif,)

    os.chdir(origdir)
    shutil.rmtree(workdir)

    result = " ".join(str(x) for x in params)
    print result

    if(args.output != ""):
        with open(os.path.join(args.output, name.replace(".root", ".txt")), "w") as f:
            f.write(result)
