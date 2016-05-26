#! /usr/bin/env python

import argparse
import functools
import os
import ROOT
import numpy

@functools.total_ordering
class Result(object):
    def __init__(self, line):
        parts = line.split()
        self.mglu = int(parts[0])
        self.mlsp = int(parts[1])
        self.lumi = float(parts[2])
        self.met = float(parts[3])
        self.njets = int(parts[4])
        self.nbm = int(parts[5])
        self.tkveto = (parts[6] == "True")
        self.limit = float(parts[7])
        self.significance = float(parts[8])

    def __str__(self):
        return "Result(mglu="+str(self.mglu)+", mlsp="+str(self.mlsp)+", lumi="+str(self.lumi)+", tkveto="+str(self.tkveto)+", nbm="+str(self.nbm)+", njets="+str(self.njets)+", met="+str(self.met)+", limit="+str(self.limit)+", significance="+str(self.significance)+")"

    def __lt__(self, other):
        return (self.mglu, self.mlsp, self.lumi, self.tkveto, self.nbm, self.njets, self.met, self.limit, self.significance) < (other.mglu, other.mlsp, other.lumi, other.tkveto, other.nbm, other.njets, other.met, other.limit, other.significance)
    
def GroupResults(results):
    results.sort()
    out = list()
    while len(results) > 0:
        plot_match = results[0]
        plot_matches = [r for r in results if r.mglu == plot_match.mglu and r.mlsp == plot_match.mlsp and r.tkveto == plot_match.tkveto and r.nbm == plot_match.nbm]

        plot = list()
        while len(plot_matches) > 0:
            line_match = plot_matches[0]
            line_matches = [p for p in plot_matches if p.njets == line_match.njets]
        
            plot.append(line_matches)
            for l in line_matches:
                plot_matches.remove(l)
                results.remove(l)

        out.append(plot)
    return out

def GetResults(dir_name):
    results = list()
    files = [f for f in os.listdir(dir_name)]
    for f in files:
        with open(os.path.join(dir_name, f), "r") as the_file:
            for line in the_file:
                results.extend([Result(line)])

    results = GroupResults(results)
    return results

def GetMax(graphs):
    return max(max([g.GetY()[i] for i in xrange(0, g.GetN())]) for g in graphs)

def GetMin(graphs):
    return min(min([g.GetY()[i] for i in xrange(0, g.GetN())]) for g in graphs)

def GetMaxSignif(results):
    best_result = None
    for met_list in results:
        for result in met_list:
            if best_result == None or result.significance > best_result.significance:
                best_result = result
    return best_result

def GetMinLimit(results):
    best_result = None
    for met_list in results:
        for result in met_list:
            if best_result == None or result.limit < best_result.limit:
                best_result = result
    return best_result

def GetColor(index):
    if index == 0:
        return ROOT.kRed
    elif index == 1:
        return ROOT.kOrange
    elif index == 2:
        return ROOT.kGreen
    elif index == 3:
        return ROOT.kCyan
    elif index == 4:
        return ROOT.kBlue
    elif index == 5:
        return ROOT.kBlack

def GetTitle(results, is_signif):
    title =  "T1tttt("+str(results[0][0].mglu)+","+str(results[0][0].mlsp)+"), L="+str(results[0][0].lumi)+" fb^{-1}, "
    if not results[0][0].tkveto:
        title += "No "
    title += "Tk Veto, N_{b}#geq "+str(results[0][0].nbm)+": "

    if is_signif:
        title += "[signif=%(signif)0.3f]"%{"signif": GetMaxSignif(results).significance}
    else:
        title += "[limit=%(limit)0.3f]"%{"limit": GetMinLimit(results).limit}

    title += ";E_{T}^{miss} [GeV];"+("Significance" if is_signif else "X-Sec Limit/Nominal")
    return title

def Format(graphs, results, is_signif):
    the_max = GetMax(graphs)
    the_min = GetMin(graphs)
    title = GetTitle(results, is_signif)

    for i in xrange(0, len(graphs)):
        graphs[i].SetLineColor(GetColor(i))
        graphs[i].SetLineWidth(4)
        graphs[i].SetLineStyle(1)
        graphs[i].SetMarkerSize(0)
        graphs[i].SetMarkerColor(0)
        graphs[i].SetFillStyle(0)
        graphs[i].SetFillColor(0)
        graphs[i].SetTitle(title)
        graphs[i].GetHistogram().SetTitleSize(0.05, "xyz")
        graphs[i].GetHistogram().SetLabelSize(0.04, "xyz")
        graphs[i].SetMinimum(the_min)
        graphs[i].SetMaximum(the_max)

def GetFileName(ex):
    name = "t1tttt_"+str(ex.mglu)+"_"+str(ex.mlsp)+"_lumi_"+str(ex.lumi)+"_tkveto_"
    if ex.tkveto:
        name += "true"
    else:
        name += "false"
    name += "_nbm_"+str(ex.nbm)
    return name

def GetMarker(results, is_signif):
    best_result = GetMaxSignif(results) if is_signif else GetMinLimit(results)
    x = best_result.met
    y = best_result.significance if is_signif else best_result.limit
    marker = ROOT.TMarker(x,y,20)
    marker.SetMarkerSize(3)
    for i in xrange(0, len(results)):
        for j in xrange(0, len(results[i])):
            if best_result == results[i][j]:
                marker.SetMarkerColor(GetColor(i))
    return marker

def MakePlot(results, out_dir):
    sig_graphs = list()
    lim_graphs = list()

    for line in results:
        met = numpy.asarray([point.met for point in line])
        sig = numpy.asarray([point.significance for point in line])
        lim = numpy.asarray([point.limit for point in line])
        sig_graphs.extend([ROOT.TGraph(len(line), met, sig)])
        lim_graphs.extend([ROOT.TGraph(len(line), met, lim)])

    Format(sig_graphs, results, True)
    Format(lim_graphs, results, False)
                
    canvas = ROOT.TCanvas()
    canvas.SetFillStyle(4000)
    canvas.SetMargin(0.10, 0.20, 0.12, 0.12)

    legend = ROOT.TLegend(0.80, 0.12, 1.0, 0.88)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    for i in xrange(0, len(sig_graphs)):
        legend.AddEntry(sig_graphs[i], "N_{jets}#geq "+str(results[i][0].njets), "l")

    sig_marker = GetMarker(results, True)
    lim_marker = GetMarker(results, False)

    file_base = GetFileName(results[0][0])

    opt = "a l"
    for g in sig_graphs:
        g.Draw(opt)
        opt = "l same"
    legend.Draw("same")
    sig_marker.Draw("same")
    canvas.SetLogy(False)
    canvas.Print(os.path.join(out_dir, file_base+"_significance.pdf"))

    opt = "a l"
    for g in lim_graphs:
        g.Draw(opt)
        opt = "l same"
    legend.Draw("same")
    lim_marker.Draw("same")
    canvas.Print(os.path.join(out_dir, file_base+"_limit_lin.pdf"))
    canvas.SetLogy(True)
    canvas.Print(os.path.join(out_dir, file_base+"_limit_log.pdf"))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Takes a precomputed set of results with different aggregate binning and makes comparison plots")
    parser.add_argument("-i", "--input", required=True, help="Directory from which to read precomputed results")
    parser.add_argument("-o", "--output", required=True, help="Directoy in which to save plots")
    args = parser.parse_args()

    results = GetResults(args.input)
    for plot in results:
        MakePlot(plot, args.output)
