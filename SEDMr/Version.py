import subprocess
import os


def ifu_drp_version():
    cwd = os.getcwd()
    home = os.environ["HOME"]
    os.chdir(home+"/kpy/SEDMr")
    ver = subprocess.check_output(["git", "describe", "--always"]).strip()
    os.chdir(cwd)
    return ver
