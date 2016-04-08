#! /usr/bin/env python

import argparse
import tempfile
import shutil
import os
import ROOT

def ModParam(filename, workspacename, parameter, valtype, value):
    in_file = ROOT.TFile(filename, "read")
    workspace = in_file.Get(workspacename)
    param = workspace.var(parameter)
    
    if valtype == "upper":
        param.setMax(value)
    elif valtype == "lower":
        param.setMin(value)
    elif valtype == "central":
        param.setVal(value)
    else:
        print "Unknown type argument \""+valtype+".\""
        
    copy_name = tempfile.mkstemp(suffix=".root")[1]
    workspace.writeToFile(copy_name)
    in_file.Close()
    shutil.copyfile(copy_name, filename)
    os.remove(copy_name)
    ROOT.SetOwnership(in_file, True)
    ROOT.SetOwnership(workspace, True)
    ROOT.SetOwnership(param, True)
    print "Set "+valtype+" value of parameter "+parameter+" in workspace "+workspacename+" in file "+filename+" to "+str(value)+"."

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Modifies a parameter in an existing workspace.")
    parser.add_argument("-f", "--file", required=True)
    parser.add_argument("-w","--workspace", default="w")
    parser.add_argument("-p","--parameter", default="r")
    parser.add_argument("-t","--type", choices=["upper", "lower", "central"], default="upper")
    parser.add_argument("value", type=float, nargs=1)
    args = parser.parse_args()

    ModParam(args.file, args.workspace, args.parameter, args.type, args.value[0])

