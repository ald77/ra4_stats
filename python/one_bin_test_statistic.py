#! /usr/bin/env python

import argparse
import math
import numpy
import ROOT
ROOT.gROOT.SetBatch(True)

def LogPoissonLike(N, mu):
    return numpy.nan_to_num(numpy.subtract(numpy.multiply(N, numpy.log(mu)), mu))

def GetQs(mu, s, b, obs):
    muhat_uncapped = numpy.divide(numpy.subtract(obs,b), s)
    muhat_low_cap = numpy.maximum(muhat_uncapped, 0.)
    muhat_high_cap = numpy.minimum(muhat_uncapped, mu)
    muhat_capped = numpy.maximum(muhat_high_cap, 0.)
    top = LogPoissonLike(obs, b+s*mu)
    bot_standard = LogPoissonLike(obs, b)
    bot_uncapped = LogPoissonLike(obs, obs)
    bot_low_cap = LogPoissonLike(obs, numpy.add(b, numpy.multiply(s, muhat_low_cap)))
    bot_high_cap = LogPoissonLike(obs, numpy.add(b, numpy.multiply(s, muhat_high_cap)))
    bot_capped = LogPoissonLike(obs, numpy.add(b, numpy.multiply(s, muhat_capped)))
    return (obs,
            numpy.multiply(-2., numpy.subtract(top, bot_standard)),
            numpy.multiply(-2., numpy.subtract(top, bot_uncapped)),
            numpy.multiply(-2., numpy.subtract(top, bot_low_cap)),
            numpy.multiply(-2., numpy.subtract(top, bot_high_cap)),
            numpy.multiply(-2., numpy.subtract(top, bot_capped)))

def Style(h, color, style): 
    h.SetLineColor(color)
    h.SetLineStyle(style)
    h.SetLineWidth(4)

def GetName(i):
    if i == 0: return "raw_num"
    elif i == 1: return "standard"
    elif i == 2: return "no_cap"
    elif i == 3: return "low_cap"
    elif i == 4: return "high_cap"
    elif i == 5: return "both_capped"
    else: return "UNKNOWN"

def GetXmax(i, background, signal):
    if i == 0:
        r = background + signal
        return math.ceil(r+5.*math.sqrt(r))
    else:
        return 2*ROOT.RooStats.NumberCountingUtils.BinomialWithTauExpZ(signal, background, 10000.)**2 + 4

def GetXmin(i, background, signal):
    if i == 0: return max(0., math.floor(background-5.*math.sqrt(background)))
    elif i == 1: return -GetXmax(i, background, signal)
    else: return 0.

def PlotQ(i, q_0_bonly, q_0_splusb, q_1_bonly, q_1_splusb, background, signal):
    nbins = 100
    xmax = GetXmax(i, background, signal)
    xmin = GetXmin(i, background, signal)
    h_0_bonly = ROOT.TH1D("q_0_bonly", ";q;Toys", nbins, xmin, xmax)
    h_0_splusb = ROOT.TH1D("q_0_splusb", ";q;Toys", nbins, xmin, xmax)
    h_1_bonly = ROOT.TH1D("q_1_bonly", ";q;Toys", nbins, xmin, xmax)
    h_1_splusb = ROOT.TH1D("q_1_splusb", ";q;Toys", nbins, xmin, xmax)

    for q in q_0_bonly: h_0_bonly.Fill(q)
    for q in q_0_splusb: h_0_splusb.Fill(q)
    for q in q_1_bonly: h_1_bonly.Fill(q)
    for q in q_1_splusb: h_1_splusb.Fill(q)

    Style(h_0_bonly, ROOT.kBlack, 1)
    Style(h_0_splusb, ROOT.kBlack, 2)
    Style(h_1_bonly, ROOT.kRed, 1)
    Style(h_1_splusb, ROOT.kRed, 2)

    h_dumb = ROOT.TH1D("q", ";q;Toys", 1, xmin, xmax)
    h_dumb.SetMinimum(1.)
    h_dumb.SetMaximum(len(q_0_bonly))
    h_dumb.SetStats(False)
    h_dumb.GetYaxis().SetTitleOffset(1.25)

    l = ROOT.TLegend(0.1, 0.915, 0.9, 0.995)
    l.SetNColumns(4)
    l.SetBorderSize(0)
    l.AddEntry(h_0_bonly, "q_{0}^{B}", "l")
    l.AddEntry(h_0_splusb, "q_{0}^{S+B}", "l")
    l.AddEntry(h_1_bonly, "q_{1}^{B}", "l")
    l.AddEntry(h_1_splusb, "q_{1}^{S+B}", "l")

    c = ROOT.TCanvas()
    c.SetMargin(0.1, 0.05, 0.1, 0.1)
    c.SetLogy()
    h_dumb.Draw("hist")
    h_0_bonly.Draw("hist same")
    h_1_bonly.Draw("hist same")
    h_0_splusb.Draw("hist same")
    h_1_splusb.Draw("hist same")
    l.Draw()

    c.Print("plots/one_bin_test_stat_q_"+GetName(i)+"_b_"+str(background)+"_s_"+str(signal)+"_toys_"+str(len(q_0_bonly))+".pdf")

def PlotQVsN(i, n, q_0, q_1, background, signal):
    h_0 = ROOT.TH1D("", ";N_{obs};q(N_{obs})", len(n)-1, n)
    h_1 = ROOT.TH1D("", ";N_{obs};q(N_{obs})", len(n)-1, n)
    for bin in range(1, len(n)+1):
        h_0.SetBinContent(bin, q_0[bin-1])
        h_1.SetBinContent(bin, q_1[bin-1])
        h_0.SetBinError(bin, 0.)
        h_1.SetBinError(bin, 0.)

    Style(h_0, ROOT.kBlack, 1)
    Style(h_1, ROOT.kRed, 1)
    h_0.SetStats(False)
    h_1.SetStats(False)

    q_min = min(h_0.GetBinContent(h_0.GetMinimumBin()),
                h_1.GetBinContent(h_1.GetMinimumBin()))
    q_max = max(h_0.GetBinContent(h_0.GetMaximumBin()),
                h_1.GetBinContent(h_1.GetMaximumBin()))
    h_0.SetMinimum(q_min)
    h_0.SetMaximum(q_max)
    h_1.SetMinimum(q_min)
    h_1.SetMaximum(q_max)

    c = ROOT.TCanvas()
    
    h_0.Draw("hist l")
    h_1.Draw("hist l same")

    l = ROOT.TLegend(0.1, 0.915, 0.9, 0.995)
    l.SetNColumns(2)
    l.SetBorderSize(0)
    l.AddEntry(h_0, "q_{0}", "l")
    l.AddEntry(h_1, "q_{1}", "l")
    l.Draw("same")
    
    c.Print("plots/one_bin_test_stat_q_vs_n_"+GetName(i)+"_b_"+str(background)+"_s_"+str(signal)+".pdf")

def ScanSigStrength(i, mu, q_bonly, q_splusb, background, signal):
    nbins = 100
    xmax = GetXmax(i, background, signal)
    xmin = GetXmin(i, background, signal)
    h_bonly = ROOT.TH1D("", ";q_{"+str(mu)+"};Toys", nbins, xmin, xmax)
    h_splusb = ROOT.TH1D("", ";q_{"+str(mu)+"};Toys", nbins, xmin, xmax)

    for q in q_bonly: h_bonly.Fill(q)
    for q in q_splusb: h_splusb.Fill(q)

    Style(h_bonly, ROOT.kBlack, 1)
    Style(h_splusb, ROOT.kBlack, 2)

    h_dumb = ROOT.TH1D("", ";q_{"+str(mu)+"};Toys", 1, xmin, xmax)
    h_dumb.SetMinimum(1.)
    h_dumb.SetMaximum(len(q_bonly))
    h_dumb.SetStats(False)
    h_dumb.GetYaxis().SetTitleOffset(1.25)
    h_dumb.SetFillStyle(0)
    h_dumb.SetFillColor(0)
    h_dumb.SetLineStyle(0)
    h_dumb.SetLineColor(0)
    h_dumb.SetLineWidth(0)
    h_dumb.SetMarkerStyle(0)
    h_dumb.SetMarkerColor(0)
    h_dumb.SetMarkerSize(0)

    q_example = GetQs(mu, signal, background, background)[i]
    clb = numpy.mean(numpy.greater_equal(q_bonly, q_example))
    clsb = numpy.mean(numpy.greater_equal(q_splusb, q_example))
    if i == 0 or i == 1:
        clb = numpy.mean(numpy.less_equal(q_bonly, q_example))
        clsb = numpy.mean(numpy.less_equal(q_splusb, q_example))
    cls = clsb/clb

    line = ROOT.TLine(q_example, 1., q_example, len(q_bonly))
    line.SetLineColor(ROOT.kGreen)
    line.SetLineStyle(1)
    line.SetLineWidth(4)

    c = ROOT.TCanvas()
    c.SetMargin(0.1, 0.05, 0.1, 0.2)
    c.SetLogy()

    h_dumb.Draw("hist")
    h_bonly.Draw("hist same")
    h_splusb.Draw("hist same")
    line.Draw("same")

    l = ROOT.TLegend(0.1, 0.83, 0.9, 0.995)
    l.SetNColumns(2)
    l.SetBorderSize(0)
    l.SetTextSize(0.05)
    l.AddEntry(h_dumb, "r="+str(round(mu,3)), "l")
    if i!=0:
        l.AddEntry(h_splusb, "q_{"+str(round(mu,3))+"}^{rS+B} (CL_{S+B}="+str(round(clsb,3))+")", "l")
    else:
        l.AddEntry(h_splusb, "N^{rS+B} (CL_{S+B}="+str(round(clsb,3))+")", "l")
    l.AddEntry(h_dumb, "CL_{S}="+str(round(cls,3)), "l")
    if i!=0:
        l.AddEntry(h_bonly, "q_{"+str(round(mu,3))+"}^{B} (CL_{B}="+str(round(clb,3))+")", "l")
    else:
        l.AddEntry(h_bonly, "N^{B} (CL_{B}="+str(round(clb,3))+")", "l")
    l.Draw()

    c.Print("plots/one_bin_test_stat_cls_scan_"+GetName(i)+"_mu_"+str(mu)+"_b_"+str(background)+"_s_"+str(signal)+"_toys_"+str(len(q_bonly))+".pdf")
    
def OneBinTestStatistic(background, signal, num_toys):
    backgrounds = numpy.random.poisson(background, num_toys)
    splusb = numpy.random.poisson(signal+background, num_toys)

    q_0_bonly = GetQs(0., signal, background, backgrounds)
    q_0_splusb = GetQs(0., signal, background, splusb)
    q_1_bonly = GetQs(1., signal, background, backgrounds)
    q_1_splusb = GetQs(1., signal, background, splusb)
    q_2_bonly = GetQs(2., signal, background, backgrounds)
    q_2_splusb = GetQs(2., signal, background, splusb)
    q_3_bonly = GetQs(3., signal, background, backgrounds)
    q_3_splusb = GetQs(3., signal, background, splusb)
    
    nmin = GetXmin(0, background, signal)
    nmax = GetXmax(0, background, signal)
    n = numpy.linspace(nmin, nmax, 100)

    q_0 = GetQs(0., signal, background, n)
    q_1 = GetQs(1., signal, background, n)

    for i in range(0, 6):
        PlotQ(i, q_0_bonly[i], q_0_splusb[i], q_1_bonly[i], q_1_splusb[i], background, signal)
        PlotQVsN(i, n, q_0[i], q_1[i], background, signal)
    for mu in [0., 0.25, 0.5, 0.654511, 0.75, 1.]:
        obs = backgrounds
        if mu == 1.:   obs = splusb
        elif mu != 0.: obs = numpy.random.poisson(background+mu*signal, num_toys)
        q_bonly = GetQs(mu, signal, background, backgrounds)
        q_splusb = GetQs(mu, signal, background, obs)
        for i in range(0, 6):
            ScanSigStrength(i, mu, q_bonly[i], q_splusb[i], background, signal)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plots distributions of test statistics for a single bin counting experiment",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("background", type=float, help="True background rate")
    parser.add_argument("signal", type=float, help="True signal rate")
    parser.add_argument("--num_toys", type=int, default=100000, help="Number of toys to throw when plotting distributions.")
    args = parser.parse_args()

    OneBinTestStatistic(args.background, args.signal, args.num_toys)
