#! /usr/bin/env python

import argparse
import subprocess
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extracts binning parameters, limits, and significance from all single (aggregate) bin files in a directory")
    parser.add_argument("-i", "--input", default="/net/cms2/cms2r0/babymaker/wspaces/agg_bins/",
                        help="Directory containing workspace files to process")
    parser.add_argument("-o", "--output", default="agg_results",
                        help="Directory in which to store results")
    args = parser.parse_args()

    files = [f for f in os.listdir(args.input)]
    for f in files:
        full_path = os.path.join(args.input, f)
        subprocess.call(["JobSubmit.csh","./run/scan_aggregate.py","-f",full_path,"-o",args.output])
