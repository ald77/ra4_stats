import os
import subprocess
import tempfile

import ROOT

def full_path(path):
    return os.path.realpath(os.path.abspath(os.path.expanduser(path)))

def cmsenv(path):
    path = full_path(path)
    cwd = full_path(os.getcwd())
    os.chdir(path)
    with tempfile.TemporaryFile() as cout, open(os.devnull) as cerr:
        subprocess.check_call(["scramv1","runtime","-sh"],stdout=cout,stderr=cerr)
        cout.seek(0)
        for byte_line in cout:
            line = byte_line.decode("utf-8").strip()
            if line[0:6] == "unset ":
                var_list = line.split()[1:]
                for v in var_list:
                    os.unsetenv(v)
            elif line[0:7] == "export ":
                vname, vval = line[7:].rstrip(";").split("=")
                os.environ[vname] = vval.strip('"')
            else:
                raise Exception("Did not understand line: {}".format(line))
    os.chdir(cwd)

class NonROOTFileError(Exception):
    def __init__(self, path):
        self.path = path
    def __str__(self):
        return self.path+" is not a ROOT file"

class ROOTOpenError(Exception):
    def __init__(self, path, mode):
        self.path = path
        self.mode = mode
    def __str__(self):
        return "Could not open "+self.path+" in "+self.mode+" mode"

class ROOTFile(object):
    def __init__(self, path, mode):
        if os.path.splitext(path)[1] != ".root":
            raise NonROOTFileError(path)
        self.path = path
        self.mode = mode
    def __enter__(self):
        self.file = ROOT.TFile(self.path, self.mode)
        if self.file.IsZombie() or not self.file.IsOpen():
            raise ROOTOpenError(self.path, self.mode)
        return self.file
    def __exit__(self, type, value, traceback):
        self.file.Close()

