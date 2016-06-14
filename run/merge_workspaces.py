#! /usr/bin/env python

import argparse
import ROOT

def ImportPdfs(wspace_out, wspace_in, suffix, shared_vars):
    getattr(wspace_out, "import")(wspace_in.allPdfs(),
                                  ROOT.RooFit.RenameAllVariablesExcept(suffix, ",".join(shared_vars)),
                                  ROOT.RooFit.RenameAllNodes(suffix))

def ImportSet(wspace_out, wspace_2016, wspace_2015, set_name, shared_vars):
    wspace_out.defineSet(set_name, "")
    ExtendSet(wspace_out, wspace_2016, set_name, shared_vars, "2016")
    ExtendSet(wspace_out, wspace_2015, set_name, shared_vars, "2015")

def ExtendSet(wspace_out, wspace_in, set_name, shared_vars, suffix):
    iterator = wspace_in.set(set_name).createIterator()
    var = iterator.Next()
    while var:
        if var.GetName() in shared_vars:
            wspace_out.extendSet(set_name, var.GetName())
        else:
            wspace_out.extendSet(set_name, var.GetName()+"_"+suffix)
        var = iterator.Next()

def AddModel(wspace_out, config_name, model_name):
    model_config = ROOT.RooStats.ModelConfig(config_name, wspace_out)
    model_config.SetPdf(wspace_out.pdf(model_name))
    model_config.SetParametersOfInterest(wspace_out.set("POI"))
    model_config.SetObservables(wspace_out.set("observables"))
    model_config.SetNuisanceParameters(wspace_out.set("nuisances"))
    model_config.SetGlobalObservables(wspace_out.set("globalObservables"))
    getattr(wspace_out, "import")(model_config)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merges 2015 and 2016 workspaces")
    parser.add_argument("--file_2016", required=True, help="File containing 2016 workspace")
    parser.add_argument("--file_2015", required=True, help="File containing 2015 workspace")
    parser.add_argument("--file_out", required=True, help="File in which to save merged workspace")
    parser.add_argument("--wspace_2016", default="w", help="Name of 2016 workspace")
    parser.add_argument("--wspace_2015", default="w", help="Name of 2015 workspace")
    parser.add_argument("--wspace_out", default="w", help="Name of output workspace")
    parser.add_argument("-s", "--shared", nargs="*", default=["r","dilep_closure"], help="Variables shared or correlated between 2015 and 2016")
    args = parser.parse_args()

    file_2016 = ROOT.TFile(args.file_2016, "read")
    file_2015 = ROOT.TFile(args.file_2015, "read")

    wspace_2016 = file_2016.Get(args.wspace_2016)
    wspace_2015 = file_2015.Get(args.wspace_2015)
    wspace_out = ROOT.RooWorkspace(args.wspace_out)

    ImportPdfs(wspace_out, wspace_2016, "2016", args.shared)
    ImportPdfs(wspace_out, wspace_2015, "2015", args.shared)

    wspace_out.factory("PROD::model_b(model_b_2016, model_b_2015)")
    wspace_out.factory("PROD::model_s(model_s_2016, model_s_2015)")

    ImportSet(wspace_out, wspace_2016, wspace_2015, "POI", args.shared)
    ImportSet(wspace_out, wspace_2016, wspace_2015, "nuisances", args.shared)
    ImportSet(wspace_out, wspace_2016, wspace_2015, "observables", args.shared)
    ImportSet(wspace_out, wspace_2016, wspace_2015, "globalObservables", args.shared)

    data_obs = ROOT.RooDataSet("data_obs", "data_obs", wspace_out.set("observables"))
    data_obs.add(wspace_out.set("observables"))
    getattr(wspace_out, "import")(data_obs, ROOT.RooCmdArg())

    AddModel(wspace_out, "ModelConfig", "model_s")
    AddModel(wspace_out, "ModelConfig_bonly", "model_b")

    wspace_out.writeToFile(args.file_out)
