#! /usr/bin/env python

from __future__ import print_function

import argparse
import os
import math

import scipy.stats

import ROOT

import utils

def get_bin_names(w):
    funcs = w.allFunctions()
    iterator = funcs.createIterator()
    func = iterator.Next()
    names = ["Systematics"]
    while func:
        name = func.GetName()
        if name.startswith("nbkg_BLK_"):
            names.append(name[9:])
        func = iterator.Next()
    iterator.Reset()
    return names

def get_process_names(w):
    funcs = w.allFunctions()
    iterator = funcs.createIterator()
    func = iterator.Next()
    names = []
    while func:
        name = func.GetName()
        if name.startswith("nmc_BLK_") and "_PRC_" in name:
            pos = name.find("_PRC_") + 5
            names.append(name[pos:])
        func = iterator.Next()
    iterator.Reset()
    return set(names)

def get_constraint_names(w):
    pdfs = w.allPdfs()
    iterator = pdfs.createIterator()
    pdf = iterator.Next()
    names = []
    while pdf:
        name = pdf.GetName()
        if name.startswith("constraint_"):
            names.append(name[11:])
        pdf = iterator.Next()
    iterator.Reset()
    return names

def find_sat_nlls(bins, procs, w):
    for b in bins:
        name = b[0]
        b[3] = 0.
        if name == "Systematics": continue
        nobs = w.var("nobs_BLK_"+name).getVal()
        if nobs > 0:
            b[3] -= nobs*math.log(nobs) - nobs - math.lgamma(nobs+1.)
        for p in procs:
            nobs = w.var("nobsmc_BLK_{}_PRC_{}".format(name, p)).getVal()
            if nobs > 0:
                b[3] -= nobs*math.log(nobs) - nobs - math.lgamma(nobs+1.)

def set_variables(w, f):
    pars = f.floatParsFinal()
    set_r = False
    for i in range(pars.getSize()):
        par = pars[i]
        var = w.var(par.GetName())
        var.removeRange()
        var.setVal(par.getVal())
        var.setError(par.getError())
        if par.GetName() == "r":
            set_r = True
    r = w.var("r")
    if not set_r:
        r.setVal(0.)
        r.setConstant(True)
    else:
        r.setConstant(False)
        
def find_bkg_nlls(bins, procs, constraints, w, fit_b):
    set_variables(w, fit_b)
    for b in bins:
        name = b[0]
        b[1] = 0.
        if name == "Systematics":
            for c in constraints:
                b[1] -= w.pdf("constraint_{}".format(c)).getLogVal(None)
            continue
        
        b[1] -= w.pdf("pdf_null_BLK_{}".format(name)).getLogVal()
        for p in procs:
            b[1] -= w.pdf("pdf_mc_BLK_{}_PRC_{}".format(name, p)).getLogVal()
            
def find_sig_nlls(bins, procs, constraints, w, fit_s):
    set_variables(w, fit_s)
    for b in bins:
        name = b[0]
        b[2] = 0.
        if name == "Systematics":
            for c in constraints:
                b[2] -= w.pdf("constraint_{}".format(c)).getLogVal(None)
            continue

        b[2] -= w.pdf("pdf_alt_BLK_{}".format(name)).getLogVal()
        for p in procs:
            b[2] -= w.pdf("pdf_mc_BLK_{}_PRC_{}".format(name, p)).getLogVal()
            

def pretty_name(name):
    if "_BIN_" in name:
        pos = name.find("_BIN_")+5
        name = name[pos:]
    
    name = name.replace("r1_","R1: ")
    name = name.replace("r2_","R2: ")
    name = name.replace("r3_","R3: ")
    name = name.replace("r4_","R4: ")
    
    name = name.replace("lowmet", "200<MET#leq 350")
    name = name.replace("medmet", "350<MET#leq 500")
    name = name.replace("highmet", "MET>500")

    name = name.replace("_lownj", ", 6#leq N_{jets}#leq 8")
    name = name.replace("_highnj", ", N_{jets}#geq 9")

    name = name.replace("_allnb", "")
    name = name.replace("_1b", ", N_{b}=1")
    name = name.replace("_2b", ", N_{b}=2")
    name = name.replace("_3b", ", N_{b}#geq 3")

    return name

def get_threshold(n_sigma, n_bins, n_dof):
    norm_cdf = scipy.stats.norm.cdf(n_sigma)
    chi2_ppf = scipy.stats.chi2.ppf(norm_cdf, n_dof)
    return chi2_ppf/n_bins

def style(h, width, style, color):
    h.SetLineWidth(width)
    h.SetLineStyle(style)
    h.SetLineColor(color)
    
def make_plot(out_dir, input_name, bins, ndof, r_sign, bkg_idx, sig_idx):
    if sig_idx == 2:
        ndof = 1
    elif sig_idx == 3:
        r_sign = 0

    chisq = 0.
    for i in range(len(bins)):
        chisq += 2.*(bins[i][bkg_idx]-bins[i][sig_idx])
    p = scipy.stats.chi2.sf(chisq,ndof)
    if r_sign < 0:
        p = 1.-0.5*p
    elif r_sign > 0:
        p = 0.5*p
    z = scipy.stats.norm.isf(p)
    
    out_name = "{}_chisq_{}_model.pdf".format(input_name,
                                              "signal" if sig_idx == 2 else "saturated")
    out_path = os.path.join(out_dir, out_name)

    h_name = "Background-Only vs {} Model (#chi^{{2}}={:.2f}, N.D.o.F.={:d}, p={:.3f}, Z={:.2f});;#chi^{{2}}"
    h_name = h_name.format("Signal" if sig_idx == 2 else "Saturated", chisq, ndof, p, z)
                           
    h = ROOT.TH1D("", h_name, len(bins), -0.5, len(bins)-0.5)
    h.SetStats(False)
    h0 = ROOT.TH1D("", h_name, len(bins), -0.5, len(bins)-0.5)

    style(h, 4, 1, ROOT.kBlack)
    style(h0, 1, 2, ROOT.kBlack)

    y0 = get_threshold(0., len(bins), ndof)
    
    for i in range(len(bins)):
        h.GetXaxis().SetBinLabel(i+1, pretty_name(bins[i][0]))
        chisq = 2.*(bins[i][bkg_idx]-bins[i][sig_idx])
        h.SetBinContent(i+1,chisq)
        h0.SetBinContent(i+1, y0)
        
    h.LabelsOption("v","x")

    c = ROOT.TCanvas()
    c.cd()
    c.SetMargin(0.1, 0.05, 0.4, 0.1)
    h.Draw("hist")
    h0.Draw("hist same")
    c.Print(out_path)

def bin_significance(out_dir, input_path, ndof):
    input_path = utils.full_path(input_path)
    out_dir = utils.full_path(out_dir)

    input_name = os.path.splitext(os.path.basename(input_path))[0]
    
    with utils.ROOTFile(input_path, "read") as f:
        w = f.Get("w")
        fit_b = f.Get("fit_b")
        fit_s = f.Get("fit_s")
        bins = [ [b, 0., 0., 0.] for b in get_bin_names(w) ]
        procs = get_process_names(w)
        constraints = get_constraint_names(w)
        r_sign = 0
        r = fit_s.floatParsFinal().find("r")
        if r.getVal() < 0.:
            r_sign = -1
        elif r.getVal() > 0.:
            r_sign = 1

        find_sat_nlls(bins, procs, w)
        find_bkg_nlls(bins, procs, constraints, w, fit_b)
        find_sig_nlls(bins, procs, constraints, w, fit_s)

        print("{:>48s}: {:>12s} {:>12s} {:>12s}".format("Bin Name", "Bkg NLL", "Sig NLL", "Sat NLL"))
        for b in bins:
            print("{:>48s}: {:>12.3f} {:>12.3f} {:>12.3f}".format(b[0],b[1],b[2],b[3]))

        make_plot(out_dir, input_name, bins, ndof, r_sign, 1, 2)
        make_plot(out_dir, input_name, bins, ndof, 0,      1, 3)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Finds approximate significance contribution from each bin",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("in_file", help="File to analyze, containing global workspace and cached fit results")
    parser.add_argument("out_dir", default=".", nargs="?", help="Directory in which to place results")
    parser.add_argument("--ndof", type=int, default=18, help="Degrees of freedom added by saturated model compared to background-only model")

    args = parser.parse_args()

    bin_significance(args.out_dir, args.in_file, args.ndof)
