ra4_stats
=========
Statistics tools for the RA4-MJ analysis

Initial setup
-------------
ra4_stats relies on having CMSSW_7_4_14 or greater with the ROOT6 version of the Higgs combine tool installed in ~/cmssw. Only CMSSW_7_4_14 is needed to produce workspaces and compute limits and significances. Setup can be done via the command

    ./run/initial_setup.sh

# Generating workspaces

## Single workspaces
For production of a single workspace, the recommended command is

    ./run/make_workspace.exe --method m1bk -u all --lumi 2.1 --use_r4

This generates a workspace with the default setup for the compressed and non-compressed FullSim T1tttt models. This script is in the process of being deprecated in favor of run/wspace_sig.exe, which should have many of the same options, but runs on a single FastSim model point and generates three workspaces with the nominal PDF and the up and down variations.

## 2D Mass Scan
Currently, the 2D scan relies on David's batch system and only runs on the SLC6 cmsX machines. You will need to copy [Manuel's node whitelist](https://github.com/manuelfs/random/blob/master/skippednodes.list.good) into your $JOBS folder and have the batch system properly setup before proceeding. Once this is done, run

    ./run/send_sig_wspaces.py

to generate workspaces for all the points in the 2D FastSim scan.

# Getting statistical results

## Limits and significance for a single workspace
To get the observed/expected limits/significance for a single workspace, simply run

    ./run/run_combine.sh my_workspace_file.root

## Fitting for signal strength and other model parameters
Use

    ./run/extract_yields.exe -f my_workspace_file.root

to obtain maximum likelihood fit results. This will perform both a signal+background and a background-only fit to the observed yields. For both fits, it produces a table of all the fitted yields, a plot with the fitted and observed yields, a plot of the kappa/lambda factors for each bin, and a diagnostic table showing the best fit value and uncertainty on every intermediate value and parameter used in the fit model.

## 2D limit scan
Once you have all the workspaces for the 2D scan, producing limit scan plots is a two step process. The first (and by far the most time-consuming) step generates a text file containing all the observed and expected limits. This step can be run either locally or using David's batch system. Once this is done, the second step uses the text file to quickly produce a plot of the results.

### Performing the scan locally on a single computer
To extract limits using only a single computer, run

    ./run/scan.sh /Directory/Containing/Workspace/Files NumberOfJobsToRunInParallel

This will produce a single text file with all the observed and expected limits.

### Performing the scan with David's batch system
To extract limits using David's batch system, run

    ./run/send_limits.py

This will produce a single text file with all the observed and expected limits.

### Making the plot from the scan results
Both the local and batch system limit scan script produces a text file with all of the necessary observed and expected limits. Once this is done, the actual plot can be made quickly by running

    ./run/limit_scan.exe -f path_to_limits_text_file.txt

In addition to making a simple version of the limit scan plot, this script produces a .root file containing the results in a format useable by the [official limit scan tool](https://github.com/CMS-SUS-XPAG/PlotsSMS).
