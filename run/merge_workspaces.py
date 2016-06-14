#! /usr/bin/env python

import argparse
import ROOT

def GetLabel(i):
    return str(2015+i)

def ExtendSet(wspace_out, wspace_in, set_name, suffix, shared_vars):
    iterator = wspace_in.set(set_name).createIterator()
    var = iterator.Next()
    while var:
        extension = "" if var.GetName() in shared_vars else "_"+suffix
        wspace_out.extendSet(set_name, var.GetName()+extension)
        var = iterator.Next()

def AddLikelihood(wspace_out, name, num_merged):
    likelihood_list = name+"_"+(","+name+"_").join([GetLabel(i) for i in xrange(0, num_merged)])
    wspace_out.factory("PROD::"+name+"("+likelihood_list+")")

def AddModel(wspace_out, config_name, model_name):
    model_config = ROOT.RooStats.ModelConfig(config_name, wspace_out)
    model_config.SetPdf(wspace_out.pdf(model_name))
    model_config.SetParametersOfInterest(wspace_out.set("POI"))
    model_config.SetObservables(wspace_out.set("observables"))
    model_config.SetNuisanceParameters(wspace_out.set("nuisances"))
    model_config.SetGlobalObservables(wspace_out.set("globalObservables"))
    getattr(wspace_out, "import")(model_config)

def MergeWorkspaces(files_in, wspace_name_in, file_out, wspace_name_out, shared_vars):
    wspace_out = ROOT.RooWorkspace(wspace_name_out)

    var_sets = ["POI", "nuisances", "observables", "globalObservables"]
    for var_set in var_sets:
        wspace_out.defineSet(var_set, "")

    for ifile in xrange(0, len(files_in)):
        label = GetLabel(ifile)
        file_in = ROOT.TFile(files_in[ifile], "read")
        wspace_in = file_in.Get(wspace_name_in)

        getattr(wspace_out, "import")(wspace_in.allPdfs(),
                                      ROOT.RooFit.RenameAllVariablesExcept(label, ",".join(shared_vars)),
                                      ROOT.RooFit.RenameAllNodes(label))
        for var_set in var_sets:
            ExtendSet(wspace_out, wspace_in, var_set, label, shared_vars)

        file_in.Close()

    AddLikelihood(wspace_out, "model_b", len(files_in))
    AddLikelihood(wspace_out, "model_s", len(files_in))

    data_obs = ROOT.RooDataSet("data_obs", "data_obs", wspace_out.set("observables"))
    data_obs.add(wspace_out.set("observables"))
    getattr(wspace_out, "import")(data_obs, ROOT.RooCmdArg())

    AddModel(wspace_out, "ModelConfig", "model_s")
    AddModel(wspace_out, "ModelConfig_bonly", "model_b")

    wspace_out.writeToFile(file_out)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merges multiple workspaces",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("file_out", help="File in which to save merged workspace")
    parser.add_argument("file_in", nargs="+", help="File(s) containing workspace(s) to be merged. Typically 1 each for 2015 and 2016.")
    parser.add_argument("--wspace_out", default="w", help="Name of output workspace")
    parser.add_argument("--wspace_in", default="w", help="Name of input workspaces")
    parser.add_argument("-s", "--shared", nargs="*", default=["r","dilep_closure"], help="Variables shared or correlated across input workspaces")
    args = parser.parse_args()

    MergeWorkspaces(args.file_in, args.wspace_in, args.file_out, args.wspace_out, args.shared)
