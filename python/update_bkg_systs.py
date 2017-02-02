#! /usr/bin/env python

from __future__ import print_function

import argparse
import os

def fullPath(path):
    return os.path.realpath(os.path.abspath(os.path.expanduser(path)))

def ensureDir(path):
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise

def PrintBkgSysts(output_file):
    print("SYSTEMATIC dilep_lownj", file=output_file)
    print(" PROCESSES ttbar,other", file=output_file)
    print("  r2_lowmet_lownj_1b    0.06", file=output_file)
    print("  r2_lowmet_lownj_2b    0.06", file=output_file)
    print("  r2_lowmet_lownj_3b    0.06", file=output_file)
    print("  r2_medmet_lownj_1b    0.06", file=output_file)
    print("  r2_medmet_lownj_2b    0.06", file=output_file)
    print("  r2_medmet_lownj_3b    0.06", file=output_file)
    print("  r2_highmet_lownj_1b   0.06", file=output_file)
    print("  r2_highmet_lownj_2b   0.06", file=output_file)
    print("  r2_highmet_lownj_3b   0.06", file=output_file)
    print("", file=output_file)
    print("SYSTEMATIC dilep_highnj", file=output_file)
    print(" PROCESSES ttbar,other", file=output_file)
    print("  r2_lowmet_highnj_1b   0.16", file=output_file)
    print("  r2_lowmet_highnj_2b   0.16", file=output_file)
    print("  r2_lowmet_highnj_3b   0.16", file=output_file)
    print("  r2_medmet_highnj_1b   0.16", file=output_file)
    print("  r2_medmet_highnj_2b   0.16", file=output_file)
    print("  r2_medmet_highnj_3b   0.16", file=output_file)
    print("  r2_highmet_highnj_1b  0.16", file=output_file)
    print("  r2_highmet_highnj_2b  0.16", file=output_file)
    print("  r2_highmet_highnj_3b  0.16", file=output_file)
    print("", file=output_file)
    print("SYSTEMATIC fivejet_lowmet", file=output_file)
    print(" PROCESSES ttbar,other", file=output_file)
    print("  r2_lowmet_lownj_1b    0.15", file=output_file)
    print("  r2_lowmet_lownj_2b    0.15", file=output_file)
    print("  r2_lowmet_lownj_3b    0.15", file=output_file)
    print("  r2_lowmet_highnj_1b   0.15", file=output_file)
    print("  r2_lowmet_highnj_2b   0.15", file=output_file)
    print("  r2_lowmet_highnj_3b   0.15", file=output_file)
    print("", file=output_file)
    print("SYSTEMATIC fivejet_highmet", file=output_file)
    print(" PROCESSES ttbar,other", file=output_file)
    print("  r2_medmet_lownj_1b    0.37", file=output_file)
    print("  r2_medmet_lownj_2b    0.37", file=output_file)
    print("  r2_medmet_lownj_3b    0.37", file=output_file)
    print("  r2_medmet_highnj_1b   0.37", file=output_file)
    print("  r2_medmet_highnj_2b   0.37", file=output_file)
    print("  r2_medmet_highnj_3b   0.37", file=output_file)
    print("  r2_highmet_lownj_1b   0.37", file=output_file)
    print("  r2_highmet_lownj_2b   0.37", file=output_file)
    print("  r2_highmet_lownj_3b   0.37", file=output_file)
    print("  r2_highmet_highnj_1b  0.37", file=output_file)
    print("  r2_highmet_highnj_2b  0.37", file=output_file)
    print("  r2_highmet_highnj_3b  0.37", file=output_file)
    print("", file=output_file)

def UpdateOneBkgSyst(input_path, output_path, overwrite):
    print(str.join(" -> ",[input_path,output_path]))
    if os.path.exists(output_path) and not overwrite:
        raise Exception(str.join(" ", ["Output file",output_path,"already exists"]))

    with open(input_path, "r") as input_file, open(output_path, "w") as output_file:
        PrintBkgSysts(output_file)
        syst_line = ""
        proc_line = ""
        printed = False
        for line in input_file:
            if "SYSTEMATIC" in line:
                syst_line = line
                printed = False
            elif "PROCESSES" in line:
                proc_line = line
                printed = False
            elif "signal" in proc_line and syst_line != "" and proc_line != "":
                if not printed:
                    print(syst_line, file=output_file, end="")
                    print(proc_line, file=output_file, end="")
                    printed = True
                print(line, file=output_file, end="")

def UpdateBkgSysts(input_dirs, output_dir, recurse, overwrite):
    input_dirs = [ fullPath(d) for d in input_dirs ]
    output_dir = fullPath(output_dir)

    for input_dir in input_dirs:
        for root, dirs, files in os.walk(input_dir):
            if root != input_dir and not recurse:
                break
            subdir = os.path.relpath(root, input_dir)
            output_subdir = fullPath(os.path.join(output_dir, subdir))
            ensureDir(output_subdir)
            for base_file in files:
                if os.path.splitext(base_file)[1] != ".txt":
                    continue
                input_file = os.path.join(root, base_file)
                output_file = os.path.join(output_subdir, base_file)
                UpdateOneBkgSyst(input_file, output_file, overwrite)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Update systematics files in specified directories to use new background systematics.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("input_dirs", nargs="+", help="Directories containing systematics files with old background systematics")
    parser.add_argument("output_dir", help="Directory in which to place updated systematics files")
    parser.add_argument("-r","--recurse", action="store_true", help="Recurse through subdirectories")
    parser.add_argument("--overwrite", action="store_true", help="Allow overwriting of existing output files")
    args = parser.parse_args()

    UpdateBkgSysts(args.input_dirs, args.output_dir, args.recurse, args.overwrite)
