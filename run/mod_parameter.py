#! /usr/bin/env python

import argparse
import tempfile
import shutil
import os
import ROOT

parser = argparse.ArgumentParser(description="Modifies a parameter in an existing workspace.")
parser.add_argument("-f", "--file", required=True)
parser.add_argument("-w","--workspace", default="w")
parser.add_argument("-p","--parameter", default="r")
parser.add_argument("-t","--type", choices=["upper", "lower", "central"], default="upper")
parser.add_argument("value", type=float, nargs=1)
args = parser.parse_args()

in_file = ROOT.TFile(args.file, "read")
workspace = in_file.Get(args.workspace)
param = workspace.var(args.parameter)
value = args.value[0]

if args.type == "upper":
    param.setMax(value)
elif args.type == "lower":
    param.setMin(value)
elif args.type == "central":
    param.setVal(value)
else:
    print "Unknown type argument \""+args.type+".\""

copy_name = tempfile.mkstemp(suffix=".root")[1]
workspace.writeToFile(copy_name)
in_file.Close()
shutil.copyfile(copy_name, args.file)
os.remove(copy_name)
