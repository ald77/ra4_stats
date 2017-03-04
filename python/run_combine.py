#! /usr/bin/env python

from __future__ import print_function

import argparse
import os
import tempfile
import shutil
import subprocess
import sys
import math
import scipy.stats

import ROOT

import utils

ROOT.PyConfig.IgnoreCommandLineOptions = True

def run_obs_signif(rfile, overwrite, log_path):
    if not overwrite and rfile.Get("sig_obs"):
        print(" Kept observed significance: {:8.3f}".format(rfile.Get("sig_obs")[0]))
        return

    cwd = os.getcwd()
    work_dir = tempfile.mkdtemp()
    link_path = os.path.join(work_dir, os.path.basename(rfile.GetName()))
    os.symlink(rfile.GetName(), link_path)

    os.chdir(work_dir)
    command = ["combine","-M","ProfileLikelihood","--significance","--uncapped=1","--rMin=-10.","-d",link_path]
    with open(log_path, "a") as dn:
        subprocess.check_call(command, stdout=dn, stderr=dn)
    os.chdir(cwd)

    result_path = os.path.join(work_dir, "higgsCombineTest.ProfileLikelihood.mH120.root")
    signif = None
    with utils.ROOTFile(result_path, "read") as result:
        tree = result.Get("limit")
        tree.GetEntry(0)
        signif = tree.limit
        signif_vec = ROOT.TVectorD(1)
        signif_vec[0] = signif
        rfile.cd()
        signif_vec.Write("sig_obs",ROOT.TObject.kWriteDelete)
    shutil.rmtree(work_dir)
    print("Saved observed significance: {:8.3f}".format(signif))

def run_exp_signif(rfile, overwrite, log_path):
    if not overwrite and rfile.Get("sig_exp"):
        print(" Kept expected significance: {:8.3f}".format(rfile.Get("sig_exp")[0]))
        return

    cwd = os.getcwd()
    work_dir = tempfile.mkdtemp()
    link_path = os.path.join(work_dir, os.path.basename(rfile.GetName()))
    os.symlink(rfile.GetName(), link_path)

    os.chdir(work_dir)
    command = ["combine","-M","ProfileLikelihood","--significance","--uncapped=1","--rMin=-10.","--expectSignal=1","--toysFreq","-t","-1","-d",link_path]
    with open(log_path, "a") as dn:
        subprocess.check_call(command, stdout=dn, stderr=dn)
    os.chdir(cwd)

    result_path = os.path.join(work_dir, "higgsCombineTest.ProfileLikelihood.mH120.root")
    signif = None
    with utils.ROOTFile(result_path, "read") as result:
        tree = result.Get("limit")
        tree.GetEntry(0)
        signif = tree.limit
        signif_vec = ROOT.TVectorD(1)
        signif_vec[0] = signif
        rfile.cd()
        signif_vec.Write("sig_exp",ROOT.TObject.kWriteDelete)
    shutil.rmtree(work_dir)
    print("Saved expected significance: {:8.3f}".format(signif))

def run_limit(rfile, overwrite, log_path):
    if not overwrite and rfile.Get("lim_obs") and rfile.Get("lim_exp16") and rfile.Get("lim_exp50") and rfile.Get("lim_exp84"):
        print(" Kept observed limit:        {:8.3f}".format(rfile.Get("lim_obs")[0]))
        print(" Kept expected limit (16%):  {:8.3f}".format(rfile.Get("lim_exp16")[0]))
        print(" Kept expected limit (50%):  {:8.3f}".format(rfile.Get("lim_exp50")[0]))
        print(" Kept expected limit (84%):  {:8.3f}".format(rfile.Get("lim_exp84")[0]))
        return

    cwd = os.getcwd()
    work_dir = tempfile.mkdtemp()
    link_path = os.path.join(work_dir, os.path.basename(rfile.GetName()))
    os.symlink(rfile.GetName(), link_path)

    os.chdir(work_dir)
    command = ["combine","-M","Asymptotic","-d",link_path]
    with open(log_path, "a") as dn:
        subprocess.check_call(command, stdout=dn, stderr=dn)
    os.chdir(cwd)

    result_path = os.path.join(work_dir, "higgsCombineTest.Asymptotic.mH120.root")
    obs = None
    exp16 = None
    exp50 = None
    exp84 = None
    with utils.ROOTFile(result_path, "read") as result:
        tree = result.Get("limit")
        for entry in tree:
            if entry.quantileExpected < 0.:
                obs = entry.limit
            elif entry.quantileExpected > 0.15 and entry.quantileExpected < 0.17:
                exp16 = entry.limit
            elif entry.quantileExpected > 0.49 and entry.quantileExpected < 0.51:
                exp50 = entry.limit
            elif entry.quantileExpected > 0.83 and entry.quantileExpected < 0.85:
                exp84 = entry.limit
        obs_vec = ROOT.TVectorD(1)
        exp16_vec = ROOT.TVectorD(1)
        exp50_vec = ROOT.TVectorD(1)
        exp84_vec = ROOT.TVectorD(1)
        obs_vec[0] = obs
        exp16_vec[0] = exp16
        exp50_vec[0] = exp50
        exp84_vec[0] = exp84
        rfile.cd()
        obs_vec.Write("lim_obs",ROOT.TObject.kWriteDelete)
        exp16_vec.Write("lim_exp16",ROOT.TObject.kWriteDelete)
        exp50_vec.Write("lim_exp50",ROOT.TObject.kWriteDelete)
        exp84_vec.Write("lim_exp84",ROOT.TObject.kWriteDelete)
    shutil.rmtree(work_dir)
    print("Saved observed limit:        {:8.3f}".format(obs))
    print("Saved expected limit (16%):  {:8.3f}".format(exp16))
    print("Saved expected limit (50%):  {:8.3f}".format(exp50))
    print("Saved expected limit (84%):  {:8.3f}".format(exp84))

def run_fit(rfile, overwrite, log_path):
    if not overwrite and rfile.Get("fit_b") and rfile.Get("fit_s"):
        print(" Kept fit results")
        return

    cwd = os.getcwd()
    work_dir = tempfile.mkdtemp()
    link_path = os.path.join(work_dir, os.path.basename(rfile.GetName()))
    os.symlink(rfile.GetName(), link_path)

    os.chdir(work_dir)
    command = ["combine","-M","MaxLikelihoodFit","--forceRecreateNLL","--saveWorkspace","--saveWithUncertainties","--minos=all","-w","w","--dataset","data_obs","-d",link_path]
    with open(log_path, "a") as dn:
        subprocess.check_call(command, stdout=dn, stderr=dn)
    os.chdir(cwd)

    result_path = os.path.join(work_dir, "mlfit.root")
    with utils.ROOTFile(result_path, "read") as result:
        fit_b = result.Get("fit_b")
        rfile.cd()
        fit_b.Write("fit_b",ROOT.TObject.kWriteDelete)
        fit_s = result.Get("fit_s")
        rfile.cd()
        fit_s.Write("fit_s",ROOT.TObject.kWriteDelete)
    shutil.rmtree(work_dir)
    print("Saved fit results")

def fix_vars(rfile):
    if not rfile.Get("fit_b"):
        return

    w = rfile.Get("w")
    f = rfile.Get("fit_b")
    if not rfile.Get("w_orig"):
        w.Write("w_orig")

    pars = f.floatParsFinal()
    for i in range(pars.getSize()):
        p = pars.at(i)
        v = w.var(p.GetName())
        v.setVal(p.getVal())
        v.setError(p.getError())
    w.Write("w", ROOT.TObject.kWriteDelete)
    print("Set variables to background-only best fit values")

def goodness_of_fit(rfile, overwrite, ndof):
    pval = rfile.Get("pval")
    chi_sq = rfile.Get("chi2")
    saved_ndof = rfile.Get("ndof")
    if not overwrite and (pval or chi_sq or saved_ndof):
        if pval:
            print(" Kept p-value:               {:8.3f}".format(pval[0]))
        if chi_sq:
            print(" Kept chi^2:                 {:8.3f}".format(chi_sq[0]))
        if saved_ndof:
            print(" Kept N.D.o.F.:              {:4d}".format(int(saved_ndof[0])))
        return

    if not rfile.Get("w_orig"):
        print('Skipping goodness of fit. Run with "--full_fit --overwrite" first to enable.')
        return

    w = rfile.Get("w")
    pdf = w.pdf("model_b")
    fit_nll = -pdf.getLogVal()

    sat_nll = 0.
    variables = w.allVars()
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
    pval = scipy.stats.chi2.sf(chi_sq,ndof)
    chi_sq_vec = ROOT.TVectorD(1)
    ndof_vec = ROOT.TVectorD(1)
    pval_vec = ROOT.TVectorD(1)
    chi_sq_vec[0] = chi_sq
    ndof_vec[0] = float(ndof)
    pval_vec[0] = pval
    chi_sq_vec.Write("chi2",ROOT.TObject.kWriteDelete)
    ndof_vec.Write("ndof",ROOT.TObject.kWriteDelete)
    pval_vec.Write("pval",ROOT.TObject.kWriteDelete)
    
    print("Saved p-value:               {:8.3f}".format(pval))
    print("Saved chi^2:                 {:8.3f}".format(chi_sq))
    print("Saved N.D.o.F.:              {:4d}".format(int(ndof)))

def run_combine(workspace_path, output_path, do_full_fit, overwrite, log_path, ndof):
    workspace_path = utils.full_path(workspace_path)
    log_path = utils.full_path(log_path)

    with open(log_path, "w") as log_file:
        print("Workspace file: {}".format(workspace_path),file=log_file)

    if output_path is None:
        output_path = workspace_path
    else:
        output_path = utils.full_path(output_path)
        shutil.copy2(workspace_path, output_path)

    with utils.ROOTFile(output_path, "update") as rfile:
        utils.cmsenv("~/cmssw/CMSSW_7_4_14/src")
        if do_full_fit:
            run_fit(rfile, overwrite, log_path)
            if overwrite:
                fix_vars(rfile)
        goodness_of_fit(rfile, overwrite, ndof)
        run_obs_signif(rfile, overwrite, log_path)
        run_exp_signif(rfile, overwrite, log_path)
        run_limit(rfile, overwrite, log_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Runs combine to find significance, limits, and best fit parameter values",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("workspace_file", help="File containing RooWorkspace")
    parser.add_argument("--full_fit", action="store_true", help="Do full (slow) maximum likelihood fit, storing results. Also sets workspace parameters to best fit values if --overwrite option used.")
    parser.add_argument("--output_file", default=None, help="Store results in specified file instead of caching in input file")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite/update already computed results")
    parser.add_argument("--log_file", default=os.devnull, help="Log file in which to write combine output (not yet working)")
    parser.add_argument("--ndof", type=int, default=18, help="Number of degrees of freedom for goodness of fit test")

    args = parser.parse_args()

    run_combine(args.workspace_file, args.output_file, args.full_fit, args.overwrite, args.log_file, args.ndof)
