#! /usr/bin/env python

import os
import re
import mod_parameter

def GluinoMass(name):
    return int(re.search("mGluino-(.*)_mLSP-", name).group(1))

def LSPMass(name):
    return int(re.search("mLSP-(.*)_xsec", name).group(1))

workdir = "mod_sigs"
files = [f for f in os.listdir(workdir) if os.path.isfile(os.path.join(workdir,f))]

for f in files:
    mglu = GluinoMass(f)
    mlsp = LSPMass(f)

    rmax = -1
    if mglu <= 1200 and mlsp <= 550:
        rmax = 5.
        if mglu <= 900 and mlsp <= 350:
            rmax = 1.5

    if rmax > 0:
        mod_parameter.ModParam(os.path.join(workdir,f), "w", "r", "upper", rmax)
