#! /usr/bin/env python

from subprocess import call

out_dir = "/net/cms2/cms2r0/babymaker/wspaces/agg_bins"
lumi = 20.

sig_pts = ((1400,1000), (1600, 1000), (1800, 200), (1000, 300))

met_cuts = (200., 250., 300., 350., 400., 450., 500., 550., 600.)
njets_cuts = (5.5, 6.5, 7.5, 8.5, 9.5, 10.5)
nbm_cuts = (0.5, 1.5, 2.5, 3.5)

for sig_pt in sig_pts:
    mglu = sig_pt[0]
    mlsp = sig_pt[1]
    for met in met_cuts:
        for njets in njets_cuts:
            for nbm in nbm_cuts:
                call(["JobSubmit.csh", "./run/wrapper.sh", "./run/aggregate_bins.exe",
#                      "--out_dir", out_dir, "--lumi", str(lumi),
                      "--mglu", str(mglu), "--mlsp", str(mlsp),
                      "--met_low", str(met), "--met_high", str(-1.),
                      "--njets_low", str(njets), "--njets_high", str(-1.),
                      "--nbm_low", str(nbm), "--nbm_high", str(-1.),
                      "--do_track_veto", "0"])

met_cuts = (200., 350., 500.)
njets_cuts = (5.5, 8.5)
nbm_cuts = (0.5, 1.5, 2.5)

for sig_pt in sig_pts:
    mglu = sig_pt[0]
    mlsp = sig_pt[1]
    for imet in range(0, len(met_cuts)):
        met_lo = str(met_cuts[imet])
        met_hi = str(-1.)
        if imet+1 < len(met_cuts):
            met_hi = str(met_cuts[imet+1])
        for injets in range(0, len(njets_cuts)):
            njets_lo = str(njets_cuts[injets])
            njets_hi = str(-1.)
            if injets+1 < len(njets_cuts):
                njets_hi = str(njets_cuts[injets+1])
            for inbm in range(0, len(nbm_cuts)):
                nbm_lo = str(nbm_cuts[inbm])
                nbm_hi = str(-1.)
                if inbm+1 < len(nbm_cuts):
                    nbm_hi = str(nbm_cuts[inbm+1])
                call(["JobSubmit.csh", "./run/wrapper.sh", "./run/aggregate_bins.exe",
#                      "--out_dir", out_dir, "--lumi", str(lumi),
                      "--mglu", str(mglu), "--mlsp", str(mlsp),
                      "--met_low", met_lo, "--met_high", met_hi,
                      "--njets_low", njets_lo, "--njets_high", njets_hi,
                      "--nbm_low", nbm_lo, "--nbm_high", nbm_hi,
                      "--do_track_veto", "0"])
                call(["JobSubmit.csh", "./run/wrapper.sh", "./run/aggregate_bins.exe",
#                      "--out_dir", out_dir, "--lumi", str(lumi),
                      "--mglu", str(mglu), "--mlsp", str(mlsp),
                      "--met_low", met_lo, "--met_high", met_hi,
                      "--njets_low", njets_lo, "--njets_high", njets_hi,
                      "--nbm_low", nbm_lo, "--nbm_high", nbm_hi,
                      "--do_track_veto", "1"])
