#! /usr/bin/env python

from __future__ import print_function

import argparse
import numpy
import scipy.stats
import math
import ROOT

def poisson(N,mu):
    if mu>0.:
        return N*math.log(mu)-mu
    elif mu==0. and N==0.:
        return 0.
    else:
        raise ValueError("Bad poisson params: N={}, mu={}".format(N,mu))

def eval(x, coeffs):
    y = 0.
    for c in coeffs:
        y = y*x+c
    return y

class ABCDModel():
    def __init__(self, B, R):
        self.B = B
        self.R = R

    def fit(self,S,Shat = None):
        if Shat is None:
            return (self.B, self.R, S)
        elif Shat == 0.:
            T = self.B*(1.+self.R)**2+S
            Rhat = self.B*self.R*(1.+self.R)/(self.B*(1.+self.R)+S)
            Bhat = T/(1.+Rhat)**2
            return (Bhat, Rhat, Shat)
        else:
            c3 = S+self.B*(self.R+1.)**2
            c2 = -(S*S + 2.*self.B*self.R*S + 2.*self.B*S - Shat*S + self.B*self.B*self.R*self.R - 2.*Shat*self.B*self.R*self.R + 2.*self.B*self.B*self.R - 4.*Shat*self.B*self.R + self.B*self.B - Shat*self.B)
            c1 = -Shat*self.B*self.R*(2.*S + 2.*self.B*self.R - Shat*self.R + 2.*self.B - 2.*Shat)
            c0 = -Shat*Shat*self.B*self.B*self.R*self.R

            coeffs = (c3,c2,c1,c0)

            lower = self.B*self.R/(self.R+2.)
            upper = (self.B*(self.R+1.)+S)**2/(self.B*(self.R+1.)**2+S)
            Bhat = lower+0.5*(upper-lower)
            while upper>Bhat and Bhat>lower:
                if eval(Bhat, coeffs)>0.:
                    upper=Bhat
                else:
                    lower=Bhat
                Bhat = lower+0.5*(upper-lower)

            Rhat = (math.sqrt(Bhat*(4.*self.B*self.R*self.R+4.*self.B*self.R+Bhat))-Bhat)/(2.*Bhat)
            return (Bhat, Rhat, Shat)

    def nll(self, S, Bhat, Rhat, Shat):
        return -(poisson(self.B+S, Bhat+Shat)
                 +2.*poisson(self.B*self.R, Bhat*Rhat)
                 +poisson(self.B*self.R*self.R, Bhat*Rhat*Rhat))

    def q0(self, S):
        Bfix, Rfix, Sfix = self.fit(S,0.)
        Bfloat, Rfloat, Sfloat = self.fit(S)
        return 2.*(self.nll(S, Bfix, Rfix, Sfix)-self.nll(S, Bfloat, Rfloat, Sfloat))

    def qtilde(self, Stest):
        Bfix, Rfix, Sfix = self.fit(0.,Stest)
        Bfloat, Rfloat, Sfloat = self.fit(0.)
        if Sfloat > Stest:
            return 0.
        return 2.*(self.nll(0., Bfix, Rfix, Sfix)-self.nll(0., Bfloat, Rfloat, Sfloat))

    def clb(self):
        return 0.5

    def clsb(self, Stest):
        qtilde = max(self.qtilde(Stest), 0.)
        return scipy.stats.norm.sf(math.sqrt(qtilde))

    def cls(self, Stest):
        return self.clsb(Stest)/self.clb()

    def canDiscover(self, S):
        return self.q0(S) > 25.

    def canExclude(self, S):
        return self.cls(S) < 0.05

    def discoverySensitivity(self):
        lower = 0.
        upper = 1.
        while not self.canDiscover(upper):
            upper *= 2.
            if  upper >= 1000000000.:
                return upper

        mid = lower+0.5*(upper-lower)
        while upper>mid and mid>lower:
            if self.canDiscover(mid):
                upper = mid
            else:
                lower = mid
            mid = lower + 0.5*(upper-lower)
        return mid

    def upperLimit(self):
        lower = 0.
        upper = 1.
        while not self.canExclude(upper):
            upper *= 2.

        mid = lower+0.5*(upper-lower)
        while upper>mid and mid>lower:
            if self.canExclude(mid):
                upper = mid
            else:
                lower = mid
            mid = lower + 0.5*(upper-lower)
        return mid

def makeHist(bins_B, bins_R, line_style):
    h = ROOT.TH2D("",
                  ";Background Yield in Signal Region;Sideband-to-Signal Ratio;Signal Events for Exclusion/Discovery",
                  len(bins_B)-1, bins_B, len(bins_R)-1, bins_R)
    h.SetStats(False)
    h.SetMinimum(1.)
    h.SetLineWidth(5)
    h.SetLineStyle(line_style)
    h.SetLineColor(1)
    h.SetTitleOffset(1.05, "XYZ")
    h.SetTitleSize(0.045, "XYZ")
    return h

def getContours(the_max):
    one_max = math.floor(math.log(the_max,10.))
    three_max = math.floor(math.log(the_max/3.,10.))

    one_list = [ 10.**p for p in range(one_max+1) ]
    three_list = [ 3.*10.**p for p in range(one_max+1) ]

    full_list = one_list + three_list
    full_list.sort()

    return numpy.array(full_list)

def makePlot(range_B, range_R, num_bins_B, num_bins_R):
    ROOT.gStyle.SetPalette(ROOT.kRainBow)
    ROOT.gStyle.SetNumberContours(999)
    
    abcd = ABCDModel(1., 10.)
    bins_B = numpy.logspace(math.log(range_B[0],2.),
                            math.log(range_B[1],2.),
                            num_bins_B+1, base=2.)
    bins_R = numpy.logspace(math.log(range_R[0],2.),
                            math.log(range_R[1],2.),
                            num_bins_R+1, base=2.)
    h_disc = makeHist(bins_B, bins_R, 1)
    h_excl = makeHist(bins_B, bins_R, 3)
    h_dumb = makeHist(bins_B, bins_R, 1)
    
    for ib in range(num_bins_B):
        for ir in range(num_bins_R):
            abcd.B = math.sqrt(bins_B[ib]*bins_B[ib+1])
            abcd.R = math.sqrt(bins_R[ir]*bins_R[ir+1])
            h_disc.SetBinContent(ib+1, ir+1, abcd.discoverySensitivity())
            h_excl.SetBinContent(ib+1, ir+1, abcd.upperLimit())

    left_margin = 0.1
    right_margin = 0.15
    top_margin = 0.1
    bottom_margin = 0.1

    canvas = ROOT.TCanvas("canvas", "canvas", 800, 750)
    canvas.SetLogx()
    canvas.SetLogy()
    canvas.SetLogz()
    canvas.SetMargin(left_margin, right_margin, bottom_margin, top_margin)

    the_max = h_excl.GetBinContent(h_excl.GetMaximumBin())*1.001
    h_dumb.Fill(-1.,-1.)
    h_disc.SetMaximum(the_max)
    h_excl.SetMaximum(the_max)
    h_dumb.SetMaximum(the_max)

    contours = getContours(the_max)
    
    h_disc.SetContour(len(contours), contours)
    h_excl.SetContour(len(contours), contours)
    
    ROOT.gStyle.SetNumberContours(999)
    h_excl.Draw("cont1z")
    h_disc.Draw("cont1 same")

    legend = ROOT.TLegend(left_margin, 1.-top_margin, 1.-right_margin, 1.)
    legend.SetNColumns(2)
    legend.SetBorderSize(0)
    legend.AddEntry(h_disc, "5#sigma Discovery", "l")
    legend.AddEntry(h_excl, "95% C.L. Exclusion", "l")
    legend.Draw()
    
    canvas.Print("abcd_sensitivity_thresholds.pdf")
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Computes minimum number of signal events needed for ABCD sensistivity.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--rangeB", type=float, nargs=2, default=[0.1,1000.], help="Range of number of background events to test")
    parser.add_argument("--rangeR", type=float, nargs=2, default=[0.1,1000.], help="Range of sideband-to-signal ratios to test")
    parser.add_argument("--binsB", type=int, default=100, help="Number of bins to use along background axis")
    parser.add_argument("--binsR", type=int, default=100, help="Number of bins to use along ratio axis")
    args = parser.parse_args()

    makePlot(args.rangeB, args.rangeR, args.binsB, args.binsR)
