import glob
import os
import time


def report():
    """Generate DRP report using output sp_<object>.npy files"""

    flist = glob.glob("spec_*.txt")
    flist.sort()
    print("\nReport generated on %s" % time.strftime("%c"))
    print("\nSEDM DRP run in %s\nFound %d spec_*.txt files" %
          (os.getcwd(), len(flist)))
    print("See http://pharos.caltech.edu/data_access/ifu?date=%s\n" %
          os.getcwd().split('/')[-1])

    print("UTStart  Object                    Exptime Air    Flxcal"
          "   Allocation          Type Subtype  z           Rlap")
    for f in flist:
        # Get object name
        objname = f.split('_')[-1].split('.')[0]
        # Get time string
        tstr = ':'.join(f.split('_')[-4:-1])
        # check the ascii spectrum file for SNID data
        with open(f, "r") as sfl:
            lines = sfl.readlines()
            # check for SNID classification
            clas = [li for li in lines if "SNIDMATCHTYPE" in li]
            ctype = ""
            if len(clas) > 0:
                for cl in clas:
                    ctype += (" %s" % cl.split()[-1])
            # check for SNID classification
            clas = [li for li in lines if "SNIDMATCHSUBTYPE" in li]
            stype = ""
            if len(clas) > 0:
                for cl in clas:
                    stype += (" %s" % cl.split()[-1])
            # get redshift
            zmch = [li for li in lines if "REDSHIFT" in li]
            if len(zmch) > 0:
                zmch = ("%.4f" % float(zmch[0].split()[-1]))
            else:
                zmch = ""
            # get rlap
            rlap = [li for li in lines if "RLAP" in li]
            if len(rlap) > 0:
                rlap = ("%.2f" % float(rlap[0].split()[-1]))
            else:
                rlap = ""
            # get flux calibration
            flxcal = [li for li in lines if "FLUXCAL" in li]
            if len(flxcal) > 0:
                flxcal = flxcal[0].split()[-1]
            else:
                flxcal = "False"
            # get exposure time
            expt = [li for li in lines if "EXPTIME" in li]
            if len(expt) > 0:
                expt = ("%.1f" % float(expt[0].split()[-1]))
            else:
                expt = "-"
            # get airmass
            air = [li for li in lines if "AIRMASS" in li]
            if len(air) > 0:
                air = ("%.3f" % float(air[0].split()[-1]))
            else:
                air = "-"
            # get program ID
            prid = [li for li in lines if "P60PRID" in li]
            if len(prid) > 0:
                prid = prid[0].split()[-1]
            else:
                prid = "-"
            sfl.close()
        if ctype == "":
            if "STD" in f:
                ctype = " STD"

        print("%8s %-25s %7s %5s  %6s %12s  %12s  %6s %-9s  %6s" %
              (tstr, objname, expt, air, flxcal, prid, ctype, stype, zmch, rlap))


if __name__ == '__main__':
    report()
