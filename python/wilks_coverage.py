#! /usr/bin/env python

from __future__ import print_function

import argparse
import math
import scipy.stats

import ROOT

def lower_bound(N):
    if N <= 0.:
        return 0.
    
    target = N*math.log(N)-N-0.5
    bot = 0.
    top = N
    mid = bot+0.5*(top-bot)
    while top>mid and mid>bot:
        if N*math.log(mid)-mid >= target:
            top = mid
        else:
            bot = mid
        mid = bot+0.5*(top-bot)
    return mid

def upper_bound(N):
    if N <= 0.:
        return 1.
    target = N*math.log(N)-N-0.5
    bot = N
    top = max(20., N+5.*math.sqrt(N))
    mid = bot+0.5*(top-bot)
    while top>mid and mid>bot:
        if N*math.log(mid)-mid > target:
            bot = mid
        else:
            top = mid
        mid = bot+0.5*(top-bot)
    return mid

def get_bounds(upper):
    bounds = dict()
    for N in range(upper):
        bounds[N] = (lower_bound(N), upper_bound(N))
    return bounds

def get_coverage(mu, bounds):
    p = 0.
    for N in bounds:
        bound = bounds[N]
        if mu >= bound[0] and mu < bound[1]:
            p += scipy.stats.poisson.pmf(N, mu)
    return p

def test_wilks_coverage(num_pts, lower_lim, upper_lim):
    bounds = get_bounds(int(round(max(20, upper_lim+5.*math.sqrt(upper_lim)))))

    h = ROOT.TH1D("", "Wilk's Coverage;Poisson mean #mu;Coverage", num_pts, lower_lim, upper_lim)
    g = ROOT.TH1D("", "Wilk's Coverage;Poisson mean #mu;Coverage", num_pts, lower_lim, upper_lim)
    
    h.SetStats(False)
    h.SetTitleOffset(1.25, "Y")
    h.SetLineColor(ROOT.kRed)
    g.SetLineColor(ROOT.kBlack)
    g.SetLineStyle(2)

    ref_val = math.erf(1/math.sqrt(2.))
    for i in range(1, num_pts+1):
        mu = h.GetBinCenter(i)
        h.SetBinContent(i, get_coverage(mu, bounds))
        g.SetBinContent(i, ref_val)

    c = ROOT.TCanvas()
    h.Draw("L")
    g.Draw("L same")
    c.Print("wilks_coverage.pdf")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Test coverage of Poisson intervals from Wilk's theorem",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--num_pts", type=int, default=1000, help="Number of points at which to evaluate function")
    parser.add_argument("--low", type=float, default=0., help="Lower limit of tested range")
    parser.add_argument("--high", type=float, default=20., help="Upper limit of tested range")

    args = parser.parse_args()

    test_wilks_coverage(args.num_pts, args.low, args.high)
