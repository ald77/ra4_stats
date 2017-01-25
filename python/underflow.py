#! /usr/bin/env python

from __future__ import print_function

import argparse
import os
import datetime

def fullPath(path):
    return os.path.realpath(os.path.abspath(os.path.expanduser(path)))

def ensureDir(path):
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise

def Underflow(nbins, bin_yield):
    print("Underflowing datacard")
    print("date {}".format(datetime.datetime.now().strftime("%Y.%m.%d")))
    print("imax *  number of channels")
    print("jmax *  number of backgrounds")
    print("kmax *  number of nuisance parameters")
    print("------------")
    print("Observation",end="")
    for ibin in range(nbins):
        if ibin>0:
            print("\t\t{}".format(bin_yield),end="")
        else:
            print("\t{}".format(bin_yield),end="")
    print("\n------------")
    print("bin",end="")
    for ibin in range(nbins):
        print("\t{0}\t{0}".format(ibin+1),end="")
    print("\nprocess",end="")
    for ibin in range(nbins):
        print("\tS\tB",end="")
    print("\nprocess",end="")
    for ibin in range(nbins):
        print("\t0\t1",end="")
    print("\nrate",end="")
    for ibin in range(nbins):
        print("\t{0}\t{0}".format(bin_yield),end="")
    print("\n------------")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create datacard with large number of high-yield bins designed to cause underflow when fit in combine.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("nbins", type=int, help="Number of bins")
    parser.add_argument("bin_yield", type=int, help="Yield in each bin")
    args = parser.parse_args()

    Underflow(args.nbins, args.bin_yield)
