#! /usr/bin/env python

import argparse
import math
import numpy
import ROOT

def LogPoissonLike(N, mu):
    return numpy.subtract(numpy.multiply(N, numpy.log(mu)), mu)

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
    return (numpy.multiply(-2., numpy.subtract(top, bot_standard)),
            numpy.multiply(-2., numpy.subtract(top, bot_uncapped)),
            numpy.multiply(-2., numpy.subtract(top, bot_low_cap)),
            numpy.multiply(-2., numpy.subtract(top, bot_high_cap)),
            numpy.multiply(-2., numpy.subtract(top, bot_capped)))

def Style(h, color, style):
    h.SetLineColor(color)
    h.SetLineStyle(style)
    h.SetLineWidth(4)

def GetName(i):
    if i == 0: return "standard"
    elif i == 1: return "no_cap"
    elif i == 2: return "low_cap"
    elif i == 3: return "high_cap"
    elif i == 4: return "both_capped"
    else: return "UNKNOWN"

def PlotQ(i, q_0_bonly, q_0_splusb, q_1_bonly, q_1_splusb, background, signal):
    q_0_bonly = q_0_bonly[i]
    q_0_splusb = q_0_splusb[i]
    q_1_bonly = q_1_bonly[i]
    q_1_splusb = q_1_splusb[i]
    
    nbins = 100
    xmax = 2*ROOT.RooStats.NumberCountingUtils.BinomialWithTauExpZ(signal, background, 10000.)**2 + 4
    xmin = 0. if i != 0 else -xmax
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

    c.Print("q_"+GetName(i)+"_b_"+str(background)+"_s_"+str(signal)+"_toys_"+str(len(q_0_bonly))+".pdf")

def OneBinTestStatistic(background, signal, num_toys):
    backgrounds = numpy.random.poisson(background, num_toys)
    splusb = numpy.random.poisson(signal+background, num_toys)

    q_0_bonly = GetQs(0., signal, background, backgrounds)
    q_0_splusb = GetQs(0., signal, background, splusb)
    q_1_bonly = GetQs(1., signal, background, backgrounds)
    q_1_splusb = GetQs(1., signal, background, splusb)

    for i in range(0, 5):
        PlotQ(i, q_0_bonly, q_0_splusb, q_1_bonly, q_1_splusb, background, signal)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plots distributions of test statistics for a single bin counting experiment",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("background", type=float, help="True background rate")
    parser.add_argument("signal", type=float, help="True signal rate")
    parser.add_argument("--num_toys", type=int, default=100000, help="Number of toys to throw when plotting distributions.")
    args = parser.parse_args()

    OneBinTestStatistic(args.background, args.signal, args.num_toys)
