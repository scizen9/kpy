import subprocess
import os
import datetime


def ifu_drp_version():
    try:
        cwd = os.getcwd()
        home = os.environ["HOME"]
        os.chdir(home+"/kpy/SEDMr")
        ver = subprocess.check_output(["git", "describe", "--always"],
                                      shell=True).strip()
        os.chdir(cwd)
    except:
        ver = "%s" % datetime.datetime.now()

    return ver
