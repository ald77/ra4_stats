#! /usr/bin/env python

from __future__ import print_function

import argparse
import os
import sys
import math

import scipy.stats

import ROOT

def fullPath(path):
    return os.path.realpath(os.path.abspath(os.path.expanduser(path)))

class NonROOTFileError(Exception):
    def __init__(self, path):
        self.path = path
    def __str__(self):
        return self.path+" is not a ROOT file"

class ROOTOpenError(Exception):
    def __init__(self, path, mode):
        self.path = path
        self.mode = mode
    def __str__(self):
        return "Could not open "+self.path+" in "+self.mode+" mode"

class ROOTFile(object):
    def __init__(self, path, mode):
        if os.path.splitext(path)[1] != ".root":
            raise NonROOTFileError(path)
        self.path = path
        self.mode = mode
    def __enter__(self):
        self.file = ROOT.TFile(self.path, self.mode)
        if self.file.IsZombie() or not self.file.IsOpen():
            raise ROOTOpenError(self.path, self.mode)
        return self.file
    def __exit__(self, type, value, traceback):
        self.file.Close()

def GoodnessOfFit(file_path, workspace_name):
    file_path = fullPath(file_path)
    with ROOTFile(file_path, "read") as root_file:
        workspace = root_file.Get(workspace_name)
        data_obs = workspace.data("data_obs")
        model_b = workspace.pdf("model_b")

        nll = ROOT.RooNLLVar("nll", "nll", model_b, data_obs)
        nll.Print()
        minuit = ROOT.RooMinuit(nll)

        minuit.setPrintLevel(-999999)
        minuit.setVerbose(False)
        minuit.setStrategy(2)
        minuit.optimizeConst(True)
        minuit.hesse()
        minuit.migrad()
        minuit.hesse()
        minuit.migrad()
        minuit.hesse()

        fit_nll = minuit.save().minNll()

        sat_nll = 0.
        variables = workspace.allVars()
        iterator = variables.createIterator()
        variable = iterator.Next()
        while variable:
            name = variable.GetName()
            if "nobs_BLK_" in name or "nobsmc_BLK_" in name:
                val = variable.getVal()
                if val > 0:
                    sat_nll -= val*math.log(val) - val - math.lgamma(val+1.)
            variable = iterator.Next()

        chi_sq = 2.*(fit_nll-sat_nll)
        ndof = 18
        print("")
        print("ABCD NLL = {}".format(fit_nll))
        print("Saturated NLL = {}".format(sat_nll))
        print("chi^2 = {}".format(chi_sq))
        print("N.D.o.F. = {}".format(ndof))
        print("p-value = {}".format(scipy.stats.chi2.sf(chi_sq,ndof)))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Computes goodness of fit from file containing a RooWorkspace.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("input_file", help="Path to file containing RooWorkspace")
    parser.add_argument("workspace_name", default="w", nargs="?", help="Name of RooWorkspace")
    args = parser.parse_args()

    GoodnessOfFit(args.input_file, args.workspace_name)
