#!/usr/bin/env python
import os, sys, re
import glob
import string
from array import array 
from pprint import pprint
import math
import ROOT
from ROOT import TMath
from ROOT import TColor

inputfile = "master_redo.txt"

#"masteroutput5.txt"
# First, cat all output files from send_variations into a single text file
# each output file has from send_variations has a single line with an identifier, and a value (limit or significance)

def get_results_list():
    with open(inputfile) as f:
        lines = [line.rstrip('\n').split() for line in open(inputfile)] # this returns a list of pairs: [identifiying information, value of significance/limit]
        return lines


    
results = get_results_list()



#version = "mj_comparison"
version = "original_comparison"
###### Setting up plot definitions:

#Make a plot for each entry in Frames 
#Each frame has a graph for each entry defined in variations
#X axis is xvariable

#####Warning: Variations MusT NOT BE A SUBSET OF EACH OTHER: EG MJdef_mj and MJdef_mj14 . Instead use MJdef_mj_ 

xvariable = "lumi" # this will be the x axis
if "mj_comparison" in version:
    variations = [["MJdef_mj_","threshold_400"],["MJdef_mj_","threshold_350"],["MJdef_mj14","threshold_400"],["MJdef_mj14","threshold_350"]]
    var_names = ["R = 1.2, threshold = 400 GeV","R = 1.2, threshold = 350 GeV","R = 1.4, threshold = 400 GeV","R = 1.4, threshold = 350 GeV"]
    var_colors = [36,36,49,49]


if "original_comparison" in version:   
    variations = [["binning_nominal","veto_off"],["binning_nominal","veto_on"],["binning_alternate","veto_off"],["binning_alternate","veto_on"]] # Each list will define one set of limits/sigs to be plotted as one TGraph
    var_names = ["RA4 2015","RA4 2015 + veto","Alternate binning","Alternate binning + veto"] #Legend names for each TGraph
    var_colors =[30,30,46,46]

    
var_linestyles = [1,7,1,7]
var_markerstyles = [20,22,20,22]
frames = [["mGluino-1600_mLSP-1000","sig_str0"],["mGluino-1600_mLSP-1000","sig_str1"],["mGluino-1800_mLSP-200","sig_str0"],["mGluino-1800_mLSP-200","sig_str1"],["mGluino-1400_mLSP-1000","sig_str0"],["mGluino-1400_mLSP-1000","sig_str1"]] # limit, sig plots for compressed, noncompressed 
frame_titles = ["T1tttt(1600,1000) Limit","T1tttt(1600,1000) Significance","T1tttt(1800,200) Limit","T1tttt(1800,200) Significance","T1tttt(1400,1000) Limit","T1tttt(1400,1000) Significance"]

### for Full sim points:
#frames = [["mGluino-1200_mLSP-800","sig_str0"],["mGluino-1200_mLSP-800","sig_str1"],["mGluino-1500_mLSP-100","sig_str0"],["mGluino-1500_mLSP-100","sig_str1"]] # limit, sig plots for compressed, noncompressed 
#frame_titles = ["T1tttt(1200,800) Limit","T1tttt(1200,800) Significance","T1tttt(1500,100) Limit","T1tttt(1500,100) Significance"]



#####Cosmetics

ROOT.gROOT.SetBatch(ROOT.kTRUE) #prevent th1->Draw() from trying to open X11 window
ROOT.gStyle.SetCanvasDefW(600);
ROOT.gStyle.SetCanvasDefH(600);
ROOT.gStyle.SetTitleOffset(1.2,"x") 
ROOT.gStyle.SetTitleOffset(1.7,"y")
ROOT.gStyle.SetPadLeftMargin(0.14)
ROOT.gStyle.SetPadBottomMargin(0.14)
ROOT.gStyle.SetLabelFont(42)
ROOT.gStyle.SetTitleFont(42)

#######

#####Reading values and plotting them

xmin=4
xmax=21
for ifr,frame in enumerate(frames):
    graphs = []
    for var in variations:
        pairs = []
        for line in results:
            use_this_line = True
            for attribute in frame: # Find out if this line belongs in this frame
                if attribute not in line[0]:
                    use_this_line = False
                    break
            if use_this_line: #Correct frame, but does it belong in this variation?
                for attribute in var:
                    if attribute not in line[0]:
                        use_this_line = False
                        break

                if use_this_line:
                    yval = line[1]
                    xval = line[0].split(xvariable)[1].split("_")[0].split("ifb")[0] # this finds the string in the identifier between "xvariable" (lumi) and the next underscore. Removes trailing "ifb" just in case
                    pairs.append([float(xval),float(yval)])
                    #results.remove(line) #Each line is used exactly once, so remove this one from the list

        pairs.sort()
        x =  [pair[0] for pair in pairs] # convert from list of pairs to x list and y list for TGraph
        y =  [pair[1] for pair in pairs]
        graph = ROOT.TGraph(len(x),array("d",x),array("d",y)) #ROOT sucks
        graphs.append(graph) 
                    
    c = ROOT.TCanvas()
    leg = ROOT.TLegend(0.14,0.7,0.48,0.9)
    ymax = 3
    ymin = 0
    title = frame_titles[ifr]+";Luminosity [ifb]; Expected Significance;"
    if "sig_str0" in frame:
        ymax = 3
        ymin = 0
        title =  frame_titles[ifr]+";Luminosity [ifb]; Expected Limit;"
        
    histframe = ROOT.TH2F("frame","",10,xmin,xmax,10,ymin,ymax)
   
    histframe.SetStats(0)
    
    histframe.SetTitle(title)
    histframe.Draw()
    for i,graph in enumerate(graphs):
        graph.SetName(var_names[i])
        graph.SetLineColor(var_colors[i])
        graph.SetLineStyle(var_linestyles[i])
        graph.SetMarkerColor(var_colors[i])
        graph.SetMarkerStyle(var_markerstyles[i])
        
        graph.Draw("lp same")
        leg.AddEntry(graph,var_names[i],"lp")
        
    leg.Draw("same")
    c.Print(version+"_"+frame[0]+"_"+frame[1]+".pdf")
